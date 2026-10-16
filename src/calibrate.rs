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
//! Fit: the units are split into bins of `n` (equal numbers of units, no
//! allele number split between bins). In each bin the exponent `s_b` that
//! makes the expected number of observed units equal the observed number is
//! found by bisection. Bins whose exponents decrease with `n` are merged
//! (pool adjacent violators, re-fitting the merged bin), `S` is the exponent
//! of the last (best covered) bin and `h_b = s_b / S`, at the median allele
//! number of each bin.

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

/// Fit `S` and `h` (see the module documentation). `units` is sorted by `n`.
/// Returns the calibration (with `S_indel = S_snv`) and the final bins.
pub fn fit(units: &mut Vec<Unit>, n_bins: usize, n_min: u32) -> Result<(Calibration, Vec<Bin>)> {
    if units.iter().filter(|u| u.y).count() == 0 {
        bail!("none of the calibration sites is observed; cannot calibrate");
    }
    units.sort_by(|a, b| a.n.cmp(&b.n).then(a.r.partial_cmp(&b.r).expect("no NaN")));
    let n_bins = n_bins.max(1);
    // bin boundaries: equal numbers of units, moved so that no n is split
    let mut cuts = vec![0];
    for b in 1..n_bins {
        let mut c = units.len() * b / n_bins;
        while c < units.len() && c > 0 && units[c].n == units[c - 1].n {
            c += 1;
        }
        if c > *cuts.last().unwrap() && c < units.len() {
            cuts.push(c);
        }
    }
    cuts.push(units.len());
    // blocks of consecutive bins: (start, end) in units
    let mut blocks: Vec<(usize, usize)> = cuts.windows(2).map(|w| (w[0], w[1])).collect();
    let mut bins: Vec<Bin> = blocks
        .iter()
        .map(|&(s, e)| make_bin(&units[s..e]))
        .collect();
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
    let s = bins.last().unwrap().exponent;
    if !(s.is_finite() && s > 0.0) {
        bail!(
            "cannot fit the scaling factor: the best covered bin has exponent {} \
             (every or no calibration site observed)",
            s
        );
    }
    let mut knots: Vec<(u32, f64)> = Vec::new();
    for b in &bins {
        let h = (b.exponent / s).min(1.0);
        match knots.last_mut() {
            // two bins with the same median: keep the larger h (monotone)
            Some(last) if last.0 >= b.n_median => last.1 = last.1.max(h),
            _ => knots.push((b.n_median, h)),
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
                    units.push(Unit {
                        n,
                        r,
                        y: observed.keys.contains(&key),
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
        let (c, bins) = fit(&mut units, 20, 0).unwrap();
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
        // the fit reproduces the observed counts per bin
        let mut start = 0;
        for b in &bins {
            let e = expected(&units[start..start + b.units], b.exponent);
            assert!((e - b.observed as f64).abs() < 1e-6 * b.observed as f64 + 1e-9);
            start += b.units;
        }
        // monotone
        assert!(c
            .knots()
            .windows(2)
            .all(|w| w[0].1 <= w[1].1 && w[0].0 < w[1].0));
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

    #[test]
    fn decreasing_bins_are_merged() {
        // two allele numbers, the lower one observed more often: one block
        let mut units = Vec::new();
        for i in 0..2000 {
            units.push(Unit {
                n: if i < 1000 { 1000 } else { 2000 },
                r: 1e-8,
                y: if i < 1000 { i % 5 == 0 } else { i % 10 == 0 },
                class: SubstitutionClass::Transversion,
            });
        }
        let (c, bins) = fit(&mut units, 2, 0).unwrap();
        assert_eq!(bins.len(), 1);
        assert_eq!(c.knots().len(), 1);
        assert_eq!(c.knots()[0].1, 1.0);
        assert!((expected(&units, c.s_snv) - 300.0).abs() < 1e-6);
    }

    #[test]
    fn bins_do_not_split_an_allele_number() {
        let mut units: Vec<Unit> = (0..100)
            .map(|i| Unit {
                n: if i < 70 { 500 } else { 900 },
                r: 1e-7,
                // better covered sites observed more often: no merging
                y: if i < 70 { i % 3 == 0 } else { i % 2 == 0 },
                class: SubstitutionClass::Transversion,
            })
            .collect();
        // 4 bins of 25 units would split the 70 units with n = 500
        let (_, bins) = fit(&mut units, 4, 0).unwrap();
        assert_eq!(bins.len(), 2, "{:?}", bins);
        for b in &bins {
            assert!(b.n_low == b.n_high, "{:?}", b);
        }
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
        assert!(fit(&mut units, 5, 0).is_err());
    }
}
