//! `--action calibrate`: fit the scaling factor `S` and the relative effective
//! sample size `h(n)` on (near-)neutral sites, see
//! `docs/design/allele-number-scaling.md` §5.
//!
//! Every possible synonymous point mutation (a site and an alternative base,
//! counted once even if several transcripts contain it) is a *unit* with
//! per-generation rate `r`, the allele number `n` of its site, and `y = 1` if
//! exactly that variant is among the observed variants. The model is
//! `P(y = 1) = 1 − (1 − r)^(S · h(n))`.
//!
//! Fit: the units are split into bins of `n`. The allele numbers are cut on a
//! geometric grid (steps of `an_step`, e.g. 5 %), and neighbouring grid cells
//! are joined, from the highest allele number down, until a bin has at least
//! `min_observed` observed and `min_observed` unobserved units; no allele
//! number is split between bins. Rare allele numbers (e.g. sites covered only
//! by part of the data) thus get bins of their own, however few units they
//! have compared with the best covered sites. In each bin the exponent `s_b`
//! that makes the expected number of observed units equal the observed number
//! is found by bisection. `S` is the exponent of the last (best covered) bin
//! and `h_b = s_b / S`. `h` is constant within each bin (knots at the lowest
//! and highest allele number of the bin) and linear between bins.
//!
//! `h` may decrease with `n`: the allele number also reflects data quality
//! (a poorly called region has a low allele number however many samples
//! there are), not only the sample size. With `monotone`, bins whose
//! exponents decrease with `n` are merged instead (pool adjacent violators,
//! re-fitting the merged bin), and `h` is capped at 1.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};

use mutexpect::{for_each_possible_mutation, Alteration, MutationType, SeqAnnotation};
use pattern_partition_prediction::PaPaPred;
use twobit::TwoBitFile;

use crate::allele_numbers::{AlleleNumbers, NeededSites};
use crate::io::{get_reader, get_writer};
use crate::scaling::{scale_probability, Calibration};

/// Point-mutation class, for the diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SubstitutionClass {
    CpgTransition,
    OtherTransition,
    Transversion,
}

impl SubstitutionClass {
    pub fn name(self) -> &'static str {
        match self {
            SubstitutionClass::CpgTransition => "CpG_transition",
            SubstitutionClass::OtherTransition => "other_transition",
            SubstitutionClass::Transversion => "transversion",
        }
    }

    /// Class of `reference → alt` with neighbouring bases `prev` and `next`
    /// (plus strand).
    pub fn of(prev: char, reference: char, next: char, alt: char) -> Self {
        let transition = matches!(
            (reference, alt),
            ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
        );
        if !transition {
            SubstitutionClass::Transversion
        } else if (reference == 'C' && next == 'G') || (reference == 'G' && prev == 'C') {
            SubstitutionClass::CpgTransition
        } else {
            SubstitutionClass::OtherTransition
        }
    }
}

/// One possible mutation used for the calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Unit {
    /// Allele number of the site.
    pub n: u32,
    /// Per-generation rate.
    pub r: f64,
    /// Observed?
    pub y: bool,
    pub class: SubstitutionClass,
}

/// A fitted bin of allele numbers.
#[derive(Debug, Clone, PartialEq)]
pub struct Bin {
    pub n_low: u32,
    pub n_high: u32,
    pub n_median: u32,
    pub units: usize,
    pub observed: usize,
    /// Fitted exponent (`S · h` for the bin).
    pub exponent: f64,
}

/// Expected number of observed units in `units` with exponent `s`.
fn expected(units: &[Unit], s: f64) -> f64 {
    units.iter().map(|u| -(s * (-u.r).ln_1p()).exp_m1()).sum()
}

/// The exponent `s` with `expected(units, s) = observed` (bisection on a
/// log scale). 0 if nothing is observed; `f64::INFINITY` if every unit is.
fn fit_exponent(units: &[Unit]) -> f64 {
    let observed = units.iter().filter(|u| u.y).count() as f64;
    if observed == 0.0 {
        return 0.0;
    }
    if observed >= units.len() as f64 {
        return f64::INFINITY;
    }
    let (mut low, mut high) = (1e-12f64, 1.0f64);
    while expected(units, high) < observed {
        low = high;
        high *= 16.0;
        if high > 1e30 {
            return f64::INFINITY;
        }
    }
    while expected(units, low) > observed && low > 1e-300 {
        high = low;
        low /= 16.0;
    }
    for _ in 0..200 {
        let mid = (low * high).sqrt();
        if expected(units, mid) < observed {
            low = mid;
        } else {
            high = mid;
        }
        if high / low < 1.0 + 1e-13 {
            break;
        }
    }
    (low * high).sqrt()
}

fn make_bin(units: &[Unit]) -> Bin {
    Bin {
        n_low: units[0].n,
        n_high: units[units.len() - 1].n,
        n_median: units[units.len() / 2].n,
        units: units.len(),
        observed: units.iter().filter(|u| u.y).count(),
        exponent: fit_exponent(units),
    }
}

/// How `fit` bins the units by allele number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Binning {
    /// Ratio between successive allele numbers of the grid (> 1).
    pub an_step: f64,
    /// Minimum number of observed, and of unobserved, units per bin.
    pub min_observed: usize,
    /// Merge bins whose exponents decrease with `n` (h monotone, at most 1).
    pub monotone: bool,
}

impl Default for Binning {
    fn default() -> Self {
        Binning {
            an_step: 1.05,
            min_observed: 1000,
            monotone: false,
        }
    }
}

/// Bin boundaries (start, end) in `units` (sorted by `n`), see the module
/// documentation.
fn bin_blocks(units: &[Unit], binning: &Binning) -> Vec<(usize, usize)> {
    let ln_step = binning.an_step.ln();
    let cell = |n: u32| (f64::from(n.max(1)).ln() / ln_step).floor() as i64;
    let mut blocks: Vec<(usize, usize)> = Vec::new(); // from the top down
    let mut end = units.len();
    let mut start = end;
    let (mut observed, mut unobserved) = (0usize, 0usize); // in units[start..end]
    while start > 0 {
        // extend the current bin by the next lower grid cell
        let c = cell(units[start - 1].n);
        while start > 0 && cell(units[start - 1].n) == c {
            start -= 1;
            if units[start].y {
                observed += 1;
            } else {
                unobserved += 1;
            }
        }
        if observed >= binning.min_observed && unobserved >= binning.min_observed {
            blocks.push((start, end));
            end = start;
            observed = 0;
            unobserved = 0;
        }
    }
    if end > 0 {
        // the lowest units do not make a bin of their own
        match blocks.last_mut() {
            Some(last) => last.0 = 0,
            None => blocks.push((0, end)),
        }
    }
    blocks.reverse();
    blocks
}

/// Fit `S` and `h` (see the module documentation). `units` is sorted by `n`.
/// Returns the calibration (with `S_indel = S_snv`) and the final bins.
pub fn fit(units: &mut Vec<Unit>, binning: &Binning, n_min: u32) -> Result<(Calibration, Vec<Bin>)> {
    if units.iter().filter(|u| u.y).count() == 0 {
        bail!("none of the calibration sites is observed; cannot calibrate");
    }
    if !(binning.an_step > 1.0 && binning.an_step.is_finite()) {
        bail!("the allele-number step must be a number above 1, got {}", binning.an_step);
    }
    units.sort_by(|a, b| a.n.cmp(&b.n).then(a.r.partial_cmp(&b.r).expect("no NaN")));
    let mut blocks = bin_blocks(units, binning);
    let mut bins: Vec<Bin> = blocks
        .iter()
        .map(|&(s, e)| make_bin(&units[s..e]))
        .collect();
    if binning.monotone {
        // pool adjacent violators: exponents must not decrease with n
        let mut i = 0;
        while i + 1 < bins.len() {
            if bins[i].exponent > bins[i + 1].exponent {
                blocks[i].1 = blocks[i + 1].1;
                blocks.remove(i + 1);
                bins.remove(i + 1);
                bins[i] = make_bin(&units[blocks[i].0..blocks[i].1]);
                i = i.saturating_sub(1);
            } else {
                i += 1;
            }
        }
    }
    let s = bins.last().unwrap().exponent;
    if !(s.is_finite() && s > 0.0) {
        bail!(
            "cannot fit the scaling factor: the best covered bin has exponent {} \
             (every or no calibration site observed)",
            s
        );
    }
    if let Some(b) = bins.iter().find(|b| !b.exponent.is_finite()) {
        bail!(
            "every calibration site with allele number {}-{} is observed; \
             use a larger --calibration-min-observed",
            b.n_low,
            b.n_high
        );
    }
    // h is constant within a bin (knots at its lowest and highest allele
    // number), so that every bin reproduces its observed count; between bins,
    // where there are no calibration sites, it is linear
    let mut knots: Vec<(u32, f64)> = Vec::new();
    for b in &bins {
        let h = b.exponent / s;
        let h = if binning.monotone { h.min(1.0) } else { h };
        knots.push((b.n_low, h));
        if b.n_high > b.n_low {
            knots.push((b.n_high, h));
        }
    }
    let calibration = Calibration::new(s, s, n_min, knots)?;
    Ok((calibration, bins))
}

/// Observed SNVs at the needed sites, as keys `(chromosome index, 0-based
/// position, alt)`.
struct Observed {
    keys: HashSet<(u16, u32, u8)>,
}

fn read_observed(
    path: &str,
    needed: &NeededSites,
    chromosomes: &HashMap<String, u16>,
) -> Result<Observed> {
    let reader = BufReader::with_capacity(
        1 << 20,
        get_reader(path).with_context(|| format!("failed to open file {} for reading", path))?,
    );
    let mut keys = HashSet::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            bail!(
                "{}, line {}: expected at least 4 fields: chr, pos, ref, alt",
                path,
                i + 1
            );
        }
        if f[2].len() != 1 || f[3].len() != 1 {
            continue; // indel
        }
        let chr = match chromosomes.get(f[0]) {
            Some(&c) => c,
            None => continue,
        };
        let pos: usize = f[1]
            .parse()
            .with_context(|| format!("{}, line {}: bad position '{}'", path, i + 1, f[1]))?;
        if pos == 0 {
            bail!("{}, line {}: positions are 1-based", path, i + 1);
        }
        if !needed.contains(f[0], pos - 1) {
            continue;
        }
        keys.insert((
            chr,
            (pos - 1) as u32,
            f[3].as_bytes()[0].to_ascii_uppercase(),
        ));
    }
    Ok(Observed { keys })
}

/// Options of `--action calibrate`.
pub struct Options<'a> {
    pub n_min: u32,
    /// Only these transcripts (all if `None`).
    pub ids: Option<HashSet<String>>,
    pub filter_for_id: Option<&'a str>,
}

/// Collect the calibration units: possible synonymous point mutations of the
/// selected transcripts (each site and alternative base once), at sites that
/// are not masked.
pub fn collect_units(
    annotations: &[SeqAnnotation],
    genome: &TwoBitFile,
    rates: &PaPaPred,
    allele_numbers: &AlleleNumbers,
    observed_path: &str,
    needed: &NeededSites,
    options: &Options,
) -> Result<Vec<Unit>> {
    let mut chromosomes: HashMap<String, u16> = HashMap::new();
    for a in annotations {
        let next = chromosomes.len() as u16;
        chromosomes.entry(a.chr.clone()).or_insert(next);
    }
    let observed = read_observed(observed_path, needed, &chromosomes)?;
    let radius = rates.radius();
    let mut seen: HashSet<(u16, u32, u8)> = HashSet::new();
    let mut units = Vec::new();
    // per chromosome: (units, observed units), to catch observed variants
    // that do not cover all chromosomes of the annotations
    let mut per_chromosome: HashMap<String, (usize, usize)> = HashMap::new();
    let mut masked = 0usize;
    for annotation in annotations {
        if let Some(id) = options.filter_for_id {
            if annotation.name != id {
                continue;
            }
        }
        if let Some(ids) = &options.ids {
            if !ids.contains(&annotation.name) {
                continue;
            }
        }
        let start = match annotation.range.start.checked_sub(radius) {
            Some(s) => s,
            None => continue,
        };
        let seq = genome.sequence(&annotation.chr, start, annotation.range.stop + radius + 1)?;
        let chr = chromosomes[&annotation.chr];
        let result = for_each_possible_mutation(
            &seq,
            annotation,
            radius,
            false,
            false,
            false,
            |site| -> Result<()> {
                if !site
                    .mutations
                    .iter()
                    .any(|m| m.mutation_type == MutationType::Synonymous)
                {
                    return Ok(());
                }
                let probabilities = match rates.rates(site.seq_window) {
                    Ok(p) => p,
                    Err(_) => return Ok(()), // e.g. N in the context
                };
                let w: Vec<char> = site.seq_window.chars().collect();
                let n = allele_numbers.get(&annotation.chr, site.genomic_position);
                for m in site.mutations {
                    let alt = match m.alteration {
                        Alteration::Point(b) if m.mutation_type == MutationType::Synonymous => {
                            b.name()
                        }
                        _ => continue,
                    };
                    let key = (chr, site.genomic_position as u32, alt as u8);
                    if !seen.insert(key) {
                        continue; // already counted for another transcript
                    }
                    if n == 0 || n < options.n_min {
                        masked += 1;
                        continue;
                    }
                    let r = probabilities.by_name(alt) as f64;
                    if r.is_nan() {
                        continue;
                    }
                    let y = observed.keys.contains(&key);
                    let counts = per_chromosome.entry(annotation.chr.clone()).or_insert((0, 0));
                    counts.0 += 1;
                    counts.1 += y as usize;
                    units.push(Unit {
                        n,
                        r,
                        y,
                        class: SubstitutionClass::of(w[radius - 1], w[radius], w[radius + 1], alt),
                    });
                }
                Ok(())
            },
        );
        if let Err(e) = result {
            eprintln!("[WARNING] Skipping annotation {}: {}", annotation.name, e);
        }
    }
    eprintln!(
        "[INFO] calibration: {} possible synonymous mutations ({} observed), {} at masked sites left out",
        units.len(),
        units.iter().filter(|u| u.y).count(),
        masked
    );
    let mut unobserved: Vec<(&String, &(usize, usize))> = per_chromosome
        .iter()
        .filter(|(_, &(n, y))| y == 0 && n >= 1000)
        .collect();
    unobserved.sort();
    for (chr, (n, _)) in unobserved {
        eprintln!(
            "[WARNING] calibration: none of the {} possible synonymous mutations on {} is observed; \
             do the observed mutations cover {}? Otherwise leave its transcripts out (e.g. with \
             --calibration-ids or a --genomic-regions file without them)",
            n, chr, chr
        );
    }
    Ok(units)
}

/// Observed and expected counts with the fitted calibration, by allele-number
/// bin and by rate bin (deciles; separately for sites below and above the
/// median allele number), each for all units and per substitution class.
pub fn write_diagnostics(
    path: &str,
    units: &[Unit],
    bins: &[Bin],
    calibration: &Calibration,
) -> Result<()> {
    let mut out =
        get_writer(path).with_context(|| format!("failed to open file {} for writing", path))?;
    writeln!(
        out,
        "table\tbin\tan_low\tan_high\trate_low\trate_high\tclass\tunits\tobserved\texpected\tobserved_expected"
    )?;
    let p = |u: &Unit| scale_probability(u.r, calibration.s_snv * calibration.h(u.n));
    let write_rows =
        |out: &mut Box<dyn Write>, table: &str, bin: usize, group: &[&Unit]| -> Result<()> {
            if group.is_empty() {
                return Ok(());
            }
            let an_low = group.iter().map(|u| u.n).min().unwrap();
            let an_high = group.iter().map(|u| u.n).max().unwrap();
            let rate_low = group.iter().map(|u| u.r).fold(f64::INFINITY, f64::min);
            let rate_high = group.iter().map(|u| u.r).fold(0.0, f64::max);
            let classes = [
                None,
                Some(SubstitutionClass::CpgTransition),
                Some(SubstitutionClass::OtherTransition),
                Some(SubstitutionClass::Transversion),
            ];
            for class in classes.iter() {
                let members: Vec<&&Unit> = group
                    .iter()
                    .filter(|u| class.map_or(true, |c| u.class == c))
                    .collect();
                if members.is_empty() {
                    continue;
                }
                let observed = members.iter().filter(|u| u.y).count();
                let expected: f64 = members.iter().map(|u| p(u)).sum();
                writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{:e}\t{:e}\t{}\t{}\t{}\t{}\t{}",
                    table,
                    bin + 1,
                    an_low,
                    an_high,
                    rate_low,
                    rate_high,
                    class.map_or("all", |c| c.name()),
                    members.len(),
                    observed,
                    expected,
                    observed as f64 / expected
                )?;
            }
            Ok(())
        };
    // by allele number: the fitted bins (units are sorted by n)
    let mut start = 0;
    for (i, b) in bins.iter().enumerate() {
        let group: Vec<&Unit> = units[start..start + b.units].iter().collect();
        write_rows(&mut out, "by_an", i, &group)?;
        start += b.units;
    }
    // by rate, below and above the median allele number
    let median_n = units[units.len() / 2].n;
    for (table, low) in [("by_rate_low_an", true), ("by_rate_high_an", false)].iter() {
        let mut group: Vec<&Unit> = units.iter().filter(|u| (u.n < median_n) == *low).collect();
        if group.is_empty() {
            continue;
        }
        group.sort_by(|a, b| a.r.partial_cmp(&b.r).unwrap());
        for d in 0..10 {
            let slice = &group[group.len() * d / 10..group.len() * (d + 1) / 10];
            write_rows(&mut out, table, d, slice)?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Read a file with one transcript ID per line.
pub fn read_ids(path: &str) -> Result<HashSet<String>> {
    let reader = BufReader::new(
        get_reader(path).with_context(|| format!("failed to open file {} for reading", path))?,
    );
    let mut ids = HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let id = line.trim();
        if !id.is_empty() && !id.starts_with('#') {
            ids.insert(id.to_string());
        }
    }
    Ok(ids)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};

    fn true_h(n: u32) -> f64 {
        // a saturating curve with a plateau, as for exome capture kits
        let x = f64::from(n) / 100_000.0;
        if x < 1.0 {
            0.3 * x
        } else if x < 1.5 {
            0.3 + 0.1 * (x - 1.0)
        } else {
            (0.35 + 0.65 * (x - 1.5) / 1.0).min(1.0)
        }
    }

    /// Units with known S and h(n); rates spread over 3 orders of magnitude
    /// and some CpG-like high rates, so that saturation matters.
    fn simulate(s: f64, n_units: usize, seed: u64) -> Vec<Unit> {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        (0..n_units)
            .map(|_| {
                let n = rng.gen_range(5_000, 260_000);
                let cpg = rng.gen::<f64>() < 0.05;
                let r = if cpg {
                    10f64.powf(rng.gen_range(-7.3, -6.7))
                } else {
                    10f64.powf(rng.gen_range(-9.5, -7.5))
                };
                let p = scale_probability(r, s * true_h(n));
                Unit {
                    n,
                    r,
                    y: rng.gen::<f64>() < p,
                    class: if cpg {
                        SubstitutionClass::CpgTransition
                    } else {
                        SubstitutionClass::Transversion
                    },
                }
            })
            .collect()
    }

    #[test]
    fn recovers_s_and_h_from_simulated_data() {
        let s = 2.0e7;
        let mut units = simulate(s, 400_000, 1);
        let (c, bins) = fit(&mut units, &Binning::default(), 0).unwrap();
        assert!(bins.len() >= 10, "{} bins", bins.len());
        println!("S = {:e} (true {:e})", c.s_snv, s);
        for &(n, h) in c.knots() {
            println!("h({}) = {:.3} (true {:.3})", n, h, true_h(n));
        }
        // S is fitted on the best covered bin (n around 250,000, h = 1)
        assert!((c.s_snv / s - 1.0).abs() < 0.08, "S = {}", c.s_snv);
        assert_eq!(c.s_indel, c.s_snv);
        for &(n, h) in c.knots() {
            let want = true_h(n);
            assert!((h - want).abs() < 0.08, "h({}) = {}, true {}", n, h, want);
        }
        // the fit reproduces the observed counts per bin, also with h(n) of the
        // calibration (constant within a bin)
        let mut start = 0;
        for b in &bins {
            let e = expected(&units[start..start + b.units], b.exponent);
            assert!((e - b.observed as f64).abs() < 1e-6 * b.observed as f64 + 1e-9);
            let e: f64 = units[start..start + b.units]
                .iter()
                .map(|u| scale_probability(u.r, c.s_snv * c.h(u.n)))
                .sum();
            assert!((e - b.observed as f64).abs() < 1e-6 * b.observed as f64, "{:?}", b);
            start += b.units;
        }
        assert!(c.knots().windows(2).all(|w| w[0].0 < w[1].0));
        // every bin has enough observed and unobserved units
        for b in &bins {
            assert!(b.observed >= 1000 && b.units - b.observed >= 1000, "{:?}", b);
        }
    }

    #[test]
    fn exponent_matches_observed_count() {
        let units: Vec<Unit> = (0..1000)
            .map(|i| Unit {
                n: 100,
                r: 1e-8 * (1 + i % 7) as f64,
                y: i % 10 == 0,
                class: SubstitutionClass::Transversion,
            })
            .collect();
        let s = fit_exponent(&units);
        assert!((expected(&units, s) - 100.0).abs() < 1e-8);
        assert_eq!(fit_exponent(&units[1..10]), 0.0); // nothing observed
        let all: Vec<Unit> = units.iter().map(|u| Unit { y: true, ..*u }).collect();
        assert_eq!(fit_exponent(&all), f64::INFINITY);
    }

    /// Units at allele number `n` with rate 1e-8, every `every`-th observed.
    fn units_at(n: u32, count: usize, every: usize) -> Vec<Unit> {
        (0..count)
            .map(|i| Unit {
                n,
                r: 1e-8,
                y: i % every == 0,
                class: SubstitutionClass::Transversion,
            })
            .collect()
    }

    fn binning(min_observed: usize, monotone: bool) -> Binning {
        Binning {
            an_step: 1.05,
            min_observed,
            monotone,
        }
    }

    #[test]
    fn h_may_decrease_and_exceed_one() {
        // the lower allele number is observed more often (e.g. better data)
        let mut units = units_at(1000, 1000, 5);
        units.extend(units_at(2000, 1000, 10));
        let (c, bins) = fit(&mut units, &binning(50, false), 0).unwrap();
        assert_eq!(bins.len(), 2);
        assert_eq!(c.knots()[1], (2000, 1.0));
        assert!(c.knots()[0].1 > 2.0, "{:?}", c.knots());
        // each bin reproduces its own observed count
        assert!((expected(&units[..1000], c.s_snv * c.knots()[0].1) - 200.0).abs() < 1e-6);
        assert!((expected(&units[1000..], c.s_snv) - 100.0).abs() < 1e-6);
    }

    #[test]
    fn monotone_merges_decreasing_bins() {
        let mut units = units_at(1000, 1000, 5);
        units.extend(units_at(2000, 1000, 10));
        let (c, bins) = fit(&mut units, &binning(50, true), 0).unwrap();
        assert_eq!(bins.len(), 1);
        assert_eq!(c.knots(), &[(1000, 1.0), (2000, 1.0)][..]);
        assert!((expected(&units, c.s_snv) - 300.0).abs() < 1e-6);
    }

    #[test]
    fn rare_allele_numbers_get_their_own_bin() {
        // most sites fully covered, a few covered by the genomes only
        let mut units = units_at(1_600_000, 100_000, 2);
        units.extend(units_at(152_000, 2_000, 4));
        units.extend(units_at(390_000, 2_000, 3));
        let (c, bins) = fit(&mut units, &binning(200, false), 0).unwrap();
        assert_eq!(bins.len(), 3, "{:?}", bins);
        let knots: Vec<u32> = c.knots().iter().map(|k| k.0).collect();
        assert_eq!(knots, vec![152_000, 390_000, 1_600_000]);
        assert!(c.knots()[0].1 < c.knots()[1].1 && c.knots()[1].1 < 1.0);
    }

    #[test]
    fn bins_join_grid_cells_from_the_top_and_never_split_an_allele_number() {
        let mut units = units_at(1000, 300, 2); // 150 observed
        units.extend(units_at(1020, 300, 2)); // same 5% cell as 1000
        units.extend(units_at(1500, 300, 2));
        units.extend(units_at(1600, 300, 2));
        units.extend(units_at(10, 50, 2)); // too few for a bin of their own
        let (_, bins) = fit(&mut units, &binning(200, false), 0).unwrap();
        let ranges: Vec<(u32, u32)> = bins.iter().map(|b| (b.n_low, b.n_high)).collect();
        // 1600 alone has 150 observed: joined with 1500; 1000 and 1020 share a
        // cell; the 50 units at n = 10 join the lowest bin
        assert_eq!(ranges, vec![(10, 1020), (1500, 1600)], "{:?}", bins);
    }

    #[test]
    fn substitution_classes() {
        use SubstitutionClass::*;
        assert_eq!(SubstitutionClass::of('A', 'C', 'G', 'T'), CpgTransition);
        assert_eq!(SubstitutionClass::of('C', 'G', 'A', 'A'), CpgTransition);
        assert_eq!(SubstitutionClass::of('A', 'C', 'A', 'T'), OtherTransition);
        assert_eq!(SubstitutionClass::of('A', 'A', 'A', 'G'), OtherTransition);
        assert_eq!(SubstitutionClass::of('A', 'C', 'G', 'A'), Transversion);
    }

    #[test]
    fn no_observed_units_is_an_error() {
        let mut units = vec![Unit {
            n: 1,
            r: 1e-8,
            y: false,
            class: SubstitutionClass::Transversion,
        }];
        assert!(fit(&mut units, &Binning::default(), 0).is_err());
    }
}
