//! Scaling of per-generation mutation rates to the probability of observing a
//! mutation in a data set (see `docs/design/allele-number-scaling.md`, §3).
//!
//! A mutation with per-generation rate `r` at a site with allele number `n`
//! is observed with probability
//!
//! ```text
//! p = 1 − (1 − r)^(S · h(n))
//! ```
//!
//! where `S` is the scaling factor of a fully covered site (separately for
//! point mutations and indels) and `h` the relative effective sample size,
//! read from a calibration file. Sites with `n < n_min` are masked. Without a
//! calibration, `h ≡ 1` and `S = --scaling-factor`, as before.

use std::convert::TryFrom;
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};

use crate::allele_numbers::AlleleNumbers;
use crate::io::{get_reader, get_writer};
use crate::MutationType;

const FORMAT_LINE: &str = "# genovo calibration v1";

/// `1 − (1 − r)^s`, accurate also when `r` and `r · s` are tiny.
///
/// `s == 1` returns `r` unchanged; NaN stays NaN.
pub fn scale_probability(r: f64, s: f64) -> f64 {
    if s == 1.0 || r.is_nan() {
        return r;
    }
    if r >= 1.0 {
        return 1.0;
    }
    -(s * (-r).ln_1p()).exp_m1()
}

/// Which scaling factor applies to a mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    PointMutation,
    Indel,
}

impl From<MutationType> for VariantClass {
    fn from(t: MutationType) -> Self {
        match t {
            MutationType::InFrameIndel | MutationType::FrameshiftIndel => VariantClass::Indel,
            _ => VariantClass::PointMutation,
        }
    }
}

/// Scaling factors and the relative effective sample size `h(n)`.
#[derive(Debug, Clone, PartialEq)]
pub struct Calibration {
    pub s_snv: f64,
    pub s_indel: f64,
    /// Sites with a smaller allele number are masked.
    pub n_min: u32,
    /// `(n, h)`, `n` strictly increasing, `h` finite and non-negative (it may
    /// decrease and exceed 1, see `calibrate`).
    /// `h` is interpolated linearly between knots and constant outside them.
    /// No knots: `h ≡ 1`.
    knots: Vec<(u32, f64)>,
    /// Other `#` lines of the file (provenance), kept when writing.
    pub notes: Vec<String>,
}

impl Calibration {
    /// One scaling factor for everything, `h ≡ 1`, no masking: the
    /// behaviour of `--scaling-factor`.
    pub fn uniform(scaling_factor: f64) -> Self {
        Calibration {
            s_snv: scaling_factor,
            s_indel: scaling_factor,
            n_min: 0,
            knots: Vec::new(),
            notes: Vec::new(),
        }
    }

    pub fn new(s_snv: f64, s_indel: f64, n_min: u32, knots: Vec<(u32, f64)>) -> Result<Self> {
        let c = Calibration {
            s_snv,
            s_indel,
            n_min,
            knots,
            notes: Vec::new(),
        };
        c.check()?;
        Ok(c)
    }

    fn check(&self) -> Result<()> {
        for (name, s) in [("S_snv", self.s_snv), ("S_indel", self.s_indel)] {
            if !(s.is_finite() && s > 0.0) {
                bail!("{} must be a positive number, got {}", name, s);
            }
        }
        for w in self.knots.windows(2) {
            if w[1].0 <= w[0].0 {
                bail!(
                    "allele numbers of h must be increasing ({} after {})",
                    w[1].0,
                    w[0].0
                );
            }
        }
        // h may decrease with n and exceed 1 (see calibrate)
        if let Some(&(n, h)) = self.knots.iter().find(|k| !(k.1.is_finite() && k.1 >= 0.0)) {
            bail!("h must be a non-negative number, got h({}) = {}", n, h);
        }
        Ok(())
    }

    /// Whether the calibration depends on allele numbers (needs an allele
    /// number file).
    pub fn uses_allele_numbers(&self) -> bool {
        !self.knots.is_empty() || self.n_min > 0
    }

    pub fn knots(&self) -> &[(u32, f64)] {
        &self.knots
    }

    /// Whether a site with allele number `n` is masked: not callable
    /// (`n = 0`) or `n < n_min`.
    pub fn is_masked(&self, n: u32) -> bool {
        n == 0 || n < self.n_min
    }

    /// Relative effective sample size at a site with allele number `n` (not
    /// taking masking into account).
    pub fn h(&self, n: u32) -> f64 {
        let k = &self.knots;
        if k.is_empty() {
            return 1.0;
        }
        let i = k.partition_point(|&(kn, _)| kn <= n);
        if i == 0 {
            return k[0].1;
        }
        if i == k.len() {
            return k[i - 1].1;
        }
        let ((n0, h0), (n1, h1)) = (k[i - 1], k[i]);
        h0 + (h1 - h0) * f64::from(n - n0) / f64::from(n1 - n0)
    }

    fn s(&self, class: VariantClass) -> f64 {
        match class {
            VariantClass::PointMutation => self.s_snv,
            VariantClass::Indel => self.s_indel,
        }
    }

    /// The exponent `S · h(n)` for a mutation at a site with allele number
    /// `n` (`None`: no allele numbers, `h = 1`), or `None` if the site is masked.
    pub fn exponent(&self, class: VariantClass, n: Option<u32>) -> Option<f64> {
        match n {
            Some(n) if self.is_masked(n) => None,
            Some(n) => Some(self.s(class) * self.h(n)),
            None => Some(self.s(class)),
        }
    }

    pub fn read_from_file(path: &str) -> Result<Self> {
        let reader = get_reader(path)
            .with_context(|| format!("failed to open file {} for reading", path))?;
        Self::read(BufReader::new(reader)).with_context(|| format!("reading calibration {}", path))
    }

    pub fn read<R: BufRead>(reader: R) -> Result<Self> {
        let mut lines = reader.lines().enumerate();
        let first = match lines.next() {
            Some((_, line)) => line?,
            None => String::new(),
        };
        if first.trim_end() != FORMAT_LINE {
            bail!(
                "not a genovo calibration file (first line must be '{}')",
                FORMAT_LINE
            );
        }
        let (mut s_snv, mut s_indel, mut n_min) = (None, None, None);
        let mut notes = Vec::new();
        let mut knots = Vec::new();
        let mut header_seen = false;
        for (i, line) in lines {
            let line = line?;
            let text = line.trim_end();
            let line_no = i + 1;
            if text.is_empty() {
                continue;
            }
            if let Some(meta) = text.strip_prefix('#') {
                let mut fields = meta.split_whitespace();
                let key = fields.next().unwrap_or("");
                let value = fields.next();
                let parse = |v: Option<&str>| -> Result<f64> {
                    v.and_then(|v| v.parse().ok())
                        .with_context(|| format!("line {}: bad value for {}", line_no, key))
                };
                match key {
                    "S_snv" => s_snv = Some(parse(value)?),
                    "S_indel" => s_indel = Some(parse(value)?),
                    "n_min" => {
                        let v = parse(value)?;
                        n_min = Some(
                            u32::try_from(v as i64)
                                .ok()
                                .filter(|_| v.fract() == 0.0)
                                .with_context(|| {
                                    format!(
                                        "line {}: n_min must be a non-negative integer",
                                        line_no
                                    )
                                })?,
                        )
                    }
                    _ => notes.push(meta.trim().to_string()),
                }
                continue;
            }
            let fields: Vec<&str> = text.split('\t').collect();
            if !header_seen {
                if fields != ["an", "h"] {
                    bail!("line {}: expected the header 'an<TAB>h'", line_no);
                }
                header_seen = true;
                continue;
            }
            if fields.len() != 2 {
                bail!("line {}: expected 'an<TAB>h'", line_no);
            }
            let n: u32 = fields[0]
                .parse()
                .with_context(|| format!("line {}: bad allele number '{}'", line_no, fields[0]))?;
            let h: f64 = fields[1]
                .parse()
                .with_context(|| format!("line {}: bad h '{}'", line_no, fields[1]))?;
            knots.push((n, h));
        }
        let s_snv = s_snv.context("missing '# S_snv' line")?;
        let mut c = Calibration::new(s_snv, s_indel.unwrap_or(s_snv), n_min.unwrap_or(0), knots)?;
        c.notes = notes;
        Ok(c)
    }

    pub fn write<W: Write>(&self, mut out: W) -> Result<()> {
        writeln!(out, "{}", FORMAT_LINE)?;
        writeln!(out, "# S_snv\t{:e}", self.s_snv)?;
        writeln!(out, "# S_indel\t{:e}", self.s_indel)?;
        writeln!(out, "# n_min\t{}", self.n_min)?;
        for note in &self.notes {
            writeln!(out, "# {}", note)?;
        }
        writeln!(out, "an\th")?;
        for (n, h) in &self.knots {
            writeln!(out, "{}\t{}", n, h)?;
        }
        out.flush()?;
        Ok(())
    }

    pub fn write_to_file(&self, path: &str) -> Result<()> {
        let out = get_writer(path)
            .with_context(|| format!("failed to open file {} for writing", path))?;
        self.write(out)
    }
}

/// The scaling of every possible mutation: a calibration, and the allele
/// numbers if the calibration uses them.
#[derive(Debug)]
pub struct SiteScaling {
    calibration: Calibration,
    allele_numbers: Option<AlleleNumbers>,
}

impl SiteScaling {
    /// Today's behaviour: one factor for all mutations.
    #[allow(dead_code)]
    pub fn uniform(scaling_factor: f64) -> Self {
        SiteScaling {
            calibration: Calibration::uniform(scaling_factor),
            allele_numbers: None,
        }
    }

    pub fn new(calibration: Calibration, allele_numbers: Option<AlleleNumbers>) -> Result<Self> {
        if calibration.uses_allele_numbers() && allele_numbers.is_none() {
            bail!("this calibration depends on allele numbers; please provide --allele-numbers");
        }
        Ok(SiteScaling {
            calibration,
            allele_numbers,
        })
    }

    #[allow(dead_code)]
    pub fn calibration(&self) -> &Calibration {
        &self.calibration
    }

    /// Allele number at a 0-based site, if allele numbers are used.
    pub fn allele_number(&self, chr: &str, pos: usize) -> Option<u32> {
        self.allele_numbers.as_ref().map(|a| a.get(chr, pos))
    }

    /// Whether the site is masked (excluded from expected and observed counts).
    pub fn is_masked(&self, chr: &str, pos: usize) -> bool {
        self.allele_number(chr, pos)
            .map_or(false, |n| self.calibration.is_masked(n))
    }

    /// The probability of observing a mutation with per-generation rate `r` at
    /// a 0-based site, or `None` if the site is masked.
    ///
    /// For an indel, `pos` is the base after the insertion/deletion point (the
    /// position at which `enumerate` lists it; in VCF terms the anchor + 1).
    pub fn probability(&self, r: f64, class: VariantClass, chr: &str, pos: usize) -> Option<f64> {
        let n = self.allele_number(chr, pos);
        self.calibration
            .exponent(class, n)
            .map(|s| scale_probability(r, s))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele_numbers::NeededSites;
    use mutexpect::interval::Interval;
    use mutexpect::{Phase, SeqAnnotation, Strand, CDS};

    /// The formula `enumerate` used before.
    fn old_formula(r: f32, s: f32) -> f32 {
        if s == 1.0 {
            return r;
        }
        let x: f64 = 1.0 - (r as f64);
        (1.0 - x.powf(s as f64)) as f32
    }

    /// Independent reference for small r: series for ln(1 − r) and expm1.
    fn reference(r: f64, s: f64) -> f64 {
        let log1m: f64 = -(1..30).map(|k| r.powi(k) / k as f64).sum::<f64>();
        let x = s * log1m;
        if x.abs() < 1e-2 {
            // 1 − e^x = −(x + x²/2 + x³/6 + …)
            let (mut term, mut sum) = (1.0, 0.0);
            for k in 1..20 {
                term *= x / k as f64;
                sum += term;
            }
            -sum
        } else {
            1.0 - x.exp()
        }
    }

    #[test]
    fn matches_an_independent_reference() {
        for &s in &[2.0, 8586.0, 2.2e7] {
            for i in 0..2000 {
                let r = 10f64.powf(-14.0 + 11.0 * i as f64 / 2000.0); // 1e-14 .. 1e-3
                let (new, want) = (scale_probability(r, s), reference(r, s));
                assert!(
                    ((new - want) / want).abs() < 1e-12,
                    "r={} s={}: {} vs {}",
                    r,
                    s,
                    new,
                    want
                );
            }
        }
        assert_eq!(scale_probability(0.123, 1.0), 0.123); // s = 1: unchanged, as before
    }

    /// The formula `enumerate` used before (`1 − (1 − r)^s` in f64, stored as
    /// f32) agrees to f32 precision for realistic rates; it loses digits only
    /// for very small r, where `1 − r` rounds.
    #[test]
    fn close_to_the_old_formula_for_realistic_rates() {
        let old = |r: f32, s: f32| (1.0 - (1.0 - r as f64).powf(s as f64)) as f32;
        for &s in &[2.0f32, 8586.0, 2.2e7] {
            for i in 0..2000 {
                let r = 10f32.powf(-10.0 + 9.0 * i as f32 / 2000.0); // 1e-10 .. 0.1
                let new = scale_probability(r as f64, s as f64) as f32;
                let rel = ((new - old(r, s)) / new).abs();
                assert!(rel < 2e-6, "r={} s={}: {} vs {}", r, s, new, old(r, s));
            }
        }
    }

    #[test]
    fn accurate_for_tiny_probabilities() {
        // 1 − (1 − r)^s ≈ r·s for tiny r·s; the old formula loses digits here
        let (r, s) = (1e-12, 1000.0);
        let p = scale_probability(r, s);
        assert!((p - r * s).abs() / (r * s) < 1e-9, "{}", p);
        assert_eq!(scale_probability(0.0, 1e7), 0.0);
        assert_eq!(scale_probability(1.0, 1e7), 1.0);
        assert!(scale_probability(f64::NAN, 2.0).is_nan());
        assert!((scale_probability(0.5, 2.0) - 0.75).abs() < 1e-15);
    }

    fn calibration() -> Calibration {
        Calibration::new(2e7, 1e7, 1000, vec![(2000, 0.2), (4000, 0.6), (8000, 1.0)]).unwrap()
    }

    #[test]
    fn h_is_interpolated_and_constant_outside() {
        let c = calibration();
        assert_eq!(c.h(1000), 0.2);
        assert_eq!(c.h(2000), 0.2);
        assert!((c.h(3000) - 0.4).abs() < 1e-12);
        assert!((c.h(6000) - 0.8).abs() < 1e-12);
        assert_eq!(c.h(8000), 1.0);
        assert_eq!(c.h(1_000_000), 1.0);
        assert_eq!(Calibration::uniform(5.0).h(7), 1.0);
        // monotone everywhere
        let mut prev = 0.0;
        for n in (0..10000).step_by(37) {
            assert!(c.h(n) >= prev);
            prev = c.h(n);
        }
    }

    #[test]
    fn exponent_and_masking() {
        let c = calibration();
        assert_eq!(c.exponent(VariantClass::PointMutation, Some(999)), None);
        assert_eq!(
            c.exponent(VariantClass::PointMutation, Some(8000)),
            Some(2e7)
        );
        assert_eq!(c.exponent(VariantClass::Indel, Some(8000)), Some(1e7));
        assert_eq!(
            c.exponent(VariantClass::PointMutation, Some(2000)),
            Some(2e7 * 0.2)
        );
        assert_eq!(c.exponent(VariantClass::Indel, None), Some(1e7));
        assert_eq!(
            VariantClass::from(MutationType::FrameshiftIndel),
            VariantClass::Indel
        );
        assert_eq!(
            VariantClass::from(MutationType::Missense),
            VariantClass::PointMutation
        );
        assert!(c.uses_allele_numbers() && !Calibration::uniform(3.0).uses_allele_numbers());
    }

    #[test]
    fn invalid_calibrations() {
        assert!(Calibration::new(0.0, 1.0, 0, vec![]).is_err());
        assert!(Calibration::new(1.0, f64::NAN, 0, vec![]).is_err());
        assert!(Calibration::new(1.0, 1.0, 0, vec![(10, 0.5), (10, 0.6)]).is_err());
        assert!(Calibration::new(1.0, 1.0, 0, vec![(10, -0.1)]).is_err());
        assert!(Calibration::new(1.0, 1.0, 0, vec![(10, f64::NAN)]).is_err());
        // h may decrease and exceed 1
        let c = Calibration::new(1.0, 1.0, 0, vec![(10, 0.5), (20, 0.4), (30, 1.5)]).unwrap();
        assert!((c.h(15) - 0.45).abs() < 1e-12);
        assert!((c.h(25) - 0.95).abs() < 1e-12);
        assert_eq!(c.h(100), 1.5);
    }

    #[test]
    fn file_round_trip() {
        let mut c = calibration();
        c.notes = vec!["source\tsynonymous, 1000 units".to_string()];
        let mut buf = Vec::new();
        c.write(&mut buf).unwrap();
        let text = String::from_utf8(buf.clone()).unwrap();
        assert!(
            text.starts_with("# genovo calibration v1\n# S_snv\t2e7\n"),
            "{}",
            text
        );
        let back = Calibration::read(&buf[..]).unwrap();
        assert_eq!(back, c);
        // S_indel defaults to S_snv, n_min to 0, no knots: h ≡ 1
        let minimal = "# genovo calibration v1\n# S_snv 8586\nan\th\n";
        let m = Calibration::read(minimal.as_bytes()).unwrap();
        assert_eq!(m, Calibration::uniform(8586.0));
    }

    #[test]
    fn bad_files() {
        let bad = [
            "an\th\n",                                                         // no format line
            "# genovo calibration v1\nan\th\n",                                // no S_snv
            "# genovo calibration v1\n# S_snv x\nan\th\n",                     // bad value
            "# genovo calibration v1\n# S_snv 1\n100\t0.5\n",                  // no header
            "# genovo calibration v1\n# S_snv 1\nan\th\n100\t0.5\t3\n",        // extra column
            "# genovo calibration v1\n# S_snv 1\nan\th\n200\t0.5\n100\t0.6\n", // unsorted
            "# genovo calibration v1\n# S_snv 1\n# n_min 1.5\nan\th\n",        // bad n_min
        ];
        for text in bad.iter() {
            assert!(Calibration::read(text.as_bytes()).is_err(), "{}", text);
        }
    }

    fn needed_sites() -> NeededSites {
        let ann = SeqAnnotation::new(
            "T".to_string(),
            "chr1".to_string(),
            Interval::new(100, 200).unwrap(),
            Strand::Plus,
            vec![Interval::new(100, 200).unwrap()],
            vec![CDS::new(Interval::new(100, 200).unwrap(), Phase::Zero)],
        );
        NeededSites::from_annotations(&[ann], false)
    }

    #[test]
    fn site_scaling() {
        let needed = needed_sites();
        // 1-based positions 101..103: AN 500 (masked), 2000, 8000; 104 missing (0)
        let an = AlleleNumbers::read(
            "chr1\t101\t500\nchr1\t102\t2000\nchr1\t103\t8000\n".as_bytes(),
            &needed,
        )
        .unwrap();
        assert!(SiteScaling::new(calibration(), None).is_err());
        let s = SiteScaling::new(calibration(), Some(an)).unwrap();
        let r = 1e-8;
        let p = |pos| s.probability(r, VariantClass::PointMutation, "chr1", pos);
        assert_eq!(p(100), None); // masked
        assert_eq!(p(103), None); // not callable
        assert_eq!(p(101), Some(scale_probability(r, 2e7 * 0.2)));
        assert_eq!(p(102), Some(scale_probability(r, 2e7)));
        assert_eq!(
            s.probability(r, VariantClass::Indel, "chr1", 102),
            Some(scale_probability(r, 1e7))
        );
        assert!(s.is_masked("chr1", 100) && !s.is_masked("chr1", 102));
        // allele numbers without a calibration: only non-callable sites are masked
        let an = AlleleNumbers::read("chr1\t101\t500\n".as_bytes(), &needed_sites()).unwrap();
        let m = SiteScaling::new(Calibration::uniform(8586.0), Some(an)).unwrap();
        assert!(!m.is_masked("chr1", 100) && m.is_masked("chr1", 101));
        assert_eq!(
            m.probability(r, VariantClass::PointMutation, "chr1", 100),
            Some(scale_probability(r, 8586.0))
        );
        // without allele numbers: the old behaviour
        let u = SiteScaling::uniform(8586.0);
        assert_eq!(
            u.probability(r, VariantClass::Indel, "chrX", 5),
            Some(scale_probability(r, 8586.0))
        );
        assert!(!u.is_masked("chr1", 0));
    }
}
