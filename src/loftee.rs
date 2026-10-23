//! The filters of LOFTEE (Karczewski et al. 2020, <https://github.com/konradjk/loftee>,
//! `grch38` branch) for putative loss-of-function (LoF) variants: nonsense,
//! frameshift and canonical splice-site variants are high confidence (HC)
//! unless a filter marks them low confidence (LC). gnomAD v4.1 runs LOFTEE
//! without its splice predictions; of its filters, these apply to genovo's
//! possible mutations:
//!
//! * `END_TRUNC` (nonsense, frameshift): the variant is at most 50 bases
//!   upstream of the last exon-exon junction (or downstream of it), as LOFTEE
//!   measures it, and the sequence between the variant and the stop codon is
//!   poorly conserved: its GERP++ rejected-substitution scores sum to at most
//!   a cutoff (`--loftee-gerp-cutoff`; LOFTEE's default is -58, gnomAD v4.1.1
//!   uses 0).
//! * `SMALL_INTRON` (splice site): an intron shorter than 15 bases.
//! * `5UTR_SPLICE`, `3UTR_SPLICE` (splice site): an intron upstream of the
//!   start codon or downstream of the stop codon.
//!
//! Not applied: `GC_TO_GT_DONOR` (genovo's splice sites need a GT donor
//! already), `ANC_ALLELE` (needs the human ancestral sequence; rare),
//! `EXON_INTRON_UNDEF` and `INCOMPLETE_CDS` (annotation problems, the latter
//! off in gnomAD), and LOFTEE's flags, which do not make a variant LC.
//!
//! `END_TRUNC` follows LOFTEE's code (`gerp_dist.pl`, `get_last_exon_coding_length`)
//! in 1-based coordinates, including its details: the distance walked from the
//! variant to the stop codon adds up `end - start` of each exon part (not the
//! number of bases), and the exons downstream of the stop codon (3' UTR
//! exons) are added to it whole, so in a transcript with 3' UTR exons
//! `DIST_FROM_LAST_EXON` grows by their length and `END_TRUNC` rarely applies
//! (even to a variant just upstream of the last exon-exon junction). The GERP sum
//! of each part is its length times the mean score of its bases with a score,
//! as LOFTEE computes it from the bigWig file; LOFTEE takes that mean from
//! the file's summary (zoom) levels, which approximate it, so sums close to
//! the cutoff may differ.
//!
//! The coding region includes the stop codon (as the CDS of GENCODE GFF3
//! files and Ensembl's `coding_region_start/end` do). The position of a
//! frameshift is that of the base after its anchor (the first deleted base
//! of a deletion).

use std::fmt;

use mutexpect::{MutationType, SeqAnnotation, Strand};

use crate::allele_numbers::SiteValues;

/// LOFTEE's default `min_intron_size`.
pub const MIN_INTRON_SIZE: usize = 15;

/// LOFTEE's `50_BP_RULE`: `END_TRUNC` needs a distance from the last exon
/// junction of at most this.
pub const END_TRUNC_DISTANCE: i64 = 50;

/// A set of LOFTEE filters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Filters(u8);

impl Filters {
    pub const END_TRUNC: Filters = Filters(1);
    pub const SMALL_INTRON: Filters = Filters(2);
    pub const UTR5_SPLICE: Filters = Filters(4);
    pub const UTR3_SPLICE: Filters = Filters(8);

    /// Names in the order LOFTEE reports them.
    const NAMES: [(Filters, &'static str); 4] = [
        (Filters::END_TRUNC, "END_TRUNC"),
        (Filters::SMALL_INTRON, "SMALL_INTRON"),
        (Filters::UTR5_SPLICE, "5UTR_SPLICE"),
        (Filters::UTR3_SPLICE, "3UTR_SPLICE"),
    ];

    pub fn is_empty(self) -> bool {
        self.0 == 0
    }

    pub fn contains(self, other: Filters) -> bool {
        self.0 & other.0 == other.0
    }

    fn insert(&mut self, other: Filters) {
        self.0 |= other.0
    }
}

impl fmt::Display for Filters {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut first = true;
        for (filter, name) in Filters::NAMES.iter() {
            if self.contains(*filter) {
                if !first {
                    f.write_str(",")?;
                }
                f.write_str(name)?;
                first = false;
            }
        }
        Ok(())
    }
}

/// LOFTEE's verdict on one putative LoF variant, with the numbers it is
/// based on (LOFTEE's `LoF_info`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Assessment {
    pub filters: Filters,
    /// Nonsense and frameshift: GERP++ sum from the variant to the stop codon
    /// (`None` without GERP scores).
    pub gerp_dist: Option<f64>,
    /// Nonsense and frameshift: distance walked to the stop codon.
    pub bp_dist: Option<i64>,
    /// Nonsense and frameshift: `bp_dist` minus the coding length of the
    /// last coding exon (`<= 50`: the 50 bp rule fails).
    pub dist_from_last_exon: Option<i64>,
    /// Splice site: size of the intron.
    pub intron_size: Option<usize>,
}

impl Assessment {
    pub fn is_high_confidence(&self) -> bool {
        self.filters.is_empty()
    }

    /// `HC`, or `LC:` followed by the filters (as LOFTEE's `LoF` and `LoF_filter`).
    pub fn label(&self) -> String {
        if self.is_high_confidence() {
            "HC".to_string()
        } else {
            format!("LC:{}", self.filters)
        }
    }

    /// LOFTEE's `LoF_info` fields used for the filters, comma-separated
    /// (`.` if none).
    pub fn info(&self) -> String {
        let mut parts = Vec::new();
        if let Some(g) = self.gerp_dist {
            parts.push(format!("GERP_DIST:{:.3}", g));
        }
        if let Some(d) = self.bp_dist {
            parts.push(format!("BP_DIST:{}", d));
        }
        if let Some(d) = self.dist_from_last_exon {
            parts.push(format!("DIST_FROM_LAST_EXON:{}", d));
            parts.push(format!(
                "50_BP_RULE:{}",
                if d <= END_TRUNC_DISTANCE { "FAIL" } else { "PASS" }
            ));
        }
        if let Some(n) = self.intron_size {
            parts.push(format!("INTRON_SIZE:{}", n));
        }
        if parts.is_empty() {
            ".".to_string()
        } else {
            parts.join(",")
        }
    }
}

/// The settings of `--loftee`/`--loftee-filter`.
pub struct Loftee {
    gerp: Option<SiteValues<f32>>,
    gerp_cutoff: f64,
    min_intron_size: usize,
    filter: bool,
}

impl Loftee {
    /// Without `gerp` scores, `END_TRUNC` is not applied. With `filter`,
    /// low-confidence variants count as `Unknown` (see [`Loftee::counted_type`]).
    pub fn new(gerp: Option<SiteValues<f32>>, gerp_cutoff: f64, filter: bool) -> Self {
        Loftee {
            gerp,
            gerp_cutoff,
            min_intron_size: MIN_INTRON_SIZE,
            filter,
        }
    }

    pub fn filters_low_confidence(&self) -> bool {
        self.filter
    }

    /// The data needed to assess the variants of one transcript.
    pub fn transcript<'a>(&'a self, annotation: &'a SeqAnnotation) -> TranscriptLoftee<'a> {
        TranscriptLoftee::new(self, annotation)
    }

    /// The mutation type a variant of type `t` is counted as: with filtering,
    /// a low-confidence LoF variant is `Unknown` (as `--filter-plof` does with
    /// the variants it filters).
    pub fn counted_type(&self, t: MutationType, assessment: Option<&Assessment>) -> MutationType {
        match assessment {
            Some(a) if self.filter && !a.is_high_confidence() => MutationType::Unknown,
            _ => t,
        }
    }
}

/// Cumulative GERP++ sums and counts of scored bases along one exon.
struct ExonGerp {
    /// 1-based first base.
    start: i64,
    sums: Vec<f64>,
    counts: Vec<u32>,
}

impl ExonGerp {
    fn new(scores: &SiteValues<f32>, chr: &str, start: i64, end: i64) -> Self {
        let len = (end - start + 1).max(0) as usize;
        let first = (start - 1) as usize; // 0-based
        let mut sums = Vec::with_capacity(len + 1);
        let mut counts = Vec::with_capacity(len + 1);
        let (mut sum, mut count) = (0.0, 0u32);
        sums.push(sum);
        counts.push(count);
        let mut sites = scores.range(chr, first, first + len).iter().peekable();
        for k in 0..len {
            if let Some(&&(p, v)) = sites.peek() {
                if p as usize == first + k {
                    sum += f64::from(v);
                    count += 1;
                    sites.next();
                }
            }
            sums.push(sum);
            counts.push(count);
        }
        ExonGerp { start, sums, counts }
    }

    /// LOFTEE's `get_interval_gerp` of the 1-based bases `a..=b` (`a <= b`) of
    /// this exon: the number of bases times their mean score (0 without
    /// scores).
    fn interval(&self, a: i64, b: i64) -> f64 {
        let len = (self.counts.len() - 1) as i64;
        let i = (a - self.start).clamp(0, len) as usize;
        let j = (b - self.start + 1).clamp(0, len) as usize;
        if j <= i {
            return 0.0;
        }
        let count = self.counts[j] - self.counts[i];
        if count == 0 {
            return 0.0;
        }
        (b - a + 1) as f64 * (self.sums[j] - self.sums[i]) / f64::from(count)
    }
}

/// One transcript in LOFTEE's terms.
pub struct TranscriptLoftee<'a> {
    loftee: &'a Loftee,
    annotation: &'a SeqAnnotation,
    plus: bool,
    /// 1-based, inclusive, in transcript order.
    exons: Vec<(i64, i64)>,
    /// LOFTEE's `stop_codon_pos`: the 1-based coding region end on the plus
    /// strand, its start on the minus strand.
    stop: Option<i64>,
    /// 0-based, half-open range of the coding region.
    cds: Option<(usize, usize)>,
    /// Per exon (transcript order), with GERP scores.
    gerp: Option<Vec<ExonGerp>>,
}

impl<'a> TranscriptLoftee<'a> {
    fn new(loftee: &'a Loftee, annotation: &'a SeqAnnotation) -> Self {
        let plus = annotation.strand == Strand::Plus;
        let mut exons: Vec<(i64, i64)> = annotation
            .exons
            .iter()
            .map(|e| (e.start as i64 + 1, e.stop as i64))
            .collect();
        exons.sort_unstable();
        if !plus {
            exons.reverse();
        }
        let cds = annotation
            .coding_sequences
            .iter()
            .map(|c| (c.range.start, c.range.stop))
            .reduce(|(a, b), (c, d)| (a.min(c), b.max(d)));
        let stop = cds.map(|(start, stop)| if plus { stop as i64 } else { start as i64 + 1 });
        let gerp = loftee.gerp.as_ref().map(|scores| {
            exons
                .iter()
                .map(|&(s, e)| ExonGerp::new(scores, &annotation.chr, s, e))
                .collect()
        });
        TranscriptLoftee {
            loftee,
            annotation,
            plus,
            exons,
            stop,
            cds,
            gerp,
        }
    }

    /// LOFTEE's assessment of a variant of type `t` at the 0-based position
    /// `pos` (for an indel the base after its anchor); `None` if it is not a
    /// putative LoF variant.
    pub fn assess(&self, pos: usize, t: MutationType) -> Option<Assessment> {
        match t {
            MutationType::Nonsense | MutationType::FrameshiftIndel => Some(self.exonic(pos)),
            MutationType::SpliceSite => Some(self.splice(pos)),
            _ => None,
        }
    }

    /// `END_TRUNC`.
    fn exonic(&self, pos0: usize) -> Assessment {
        let mut assessment = Assessment::default();
        let stop = match self.stop {
            Some(stop) => stop,
            None => return assessment,
        };
        let pos = pos0 as i64 + 1;
        // the exon of the variant (LOFTEE: no exon number, no END_TRUNC)
        let exon_idx = match self.exons.iter().position(|&(s, e)| s <= pos && pos <= e) {
            Some(i) => i,
            None => return assessment,
        };
        // get_gerp_weighted_dist
        let mut dist = 0i64;
        let mut gerp = 0.0;
        for (k, &(s, e)) in self.exons.iter().enumerate() {
            // skip exons upstream of the variant
            if (self.plus && pos > e) || (!self.plus && pos < s) {
                continue;
            }
            let last = if self.plus {
                s < stop && e >= stop
            } else {
                e > stop && s <= stop
            };
            let affected = s <= pos && pos <= e;
            let (from, to) = if last {
                let from = if affected {
                    pos
                } else if self.plus {
                    s
                } else {
                    e
                };
                (from, stop)
            } else if affected {
                (pos, if self.plus { e } else { s })
            } else {
                (s, e)
            };
            if let Some(g) = &self.gerp {
                gerp += g[k].interval(from.min(to), from.max(to));
            }
            dist += (to - from).abs();
        }
        // get_last_exon_coding_length
        let mut last_exon_length = -1000;
        for &(s, e) in self.exons[exon_idx..].iter().rev() {
            if self.plus {
                if s > stop {
                    continue;
                } else if e >= stop {
                    last_exon_length = stop - s;
                    break;
                }
            } else if e < stop {
                continue;
            } else if s <= stop {
                last_exon_length = e - stop;
                break;
            }
        }
        let d = dist - last_exon_length;
        assessment.bp_dist = Some(dist);
        assessment.dist_from_last_exon = Some(d);
        if self.gerp.is_some() {
            assessment.gerp_dist = Some(gerp);
            if d <= END_TRUNC_DISTANCE && gerp <= self.loftee.gerp_cutoff {
                assessment.filters.insert(Filters::END_TRUNC);
            }
        }
        assessment
    }

    /// `SMALL_INTRON`, `5UTR_SPLICE`, `3UTR_SPLICE`.
    fn splice(&self, pos0: usize) -> Assessment {
        let mut assessment = Assessment::default();
        if let Some(intron) = self.annotation.find_intron(pos0) {
            assessment.intron_size = Some(intron.len());
            if intron.len() < self.loftee.min_intron_size {
                assessment.filters.insert(Filters::SMALL_INTRON);
            }
        }
        if let Some((start, stop)) = self.cds {
            let (upstream, downstream) = if self.plus {
                (pos0 < start, pos0 >= stop)
            } else {
                (pos0 >= stop, pos0 < start)
            };
            if upstream {
                assessment.filters.insert(Filters::UTR5_SPLICE);
            }
            if downstream {
                assessment.filters.insert(Filters::UTR3_SPLICE);
            }
        }
        assessment
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::allele_numbers::NeededSites;
    use mutexpect::{Interval, Phase, CDS};

    fn annotation(strand: Strand, exons: &[(usize, usize)], cds: &[(usize, usize)]) -> SeqAnnotation {
        let iv = |&(s, e): &(usize, usize)| Interval::new(s, e).unwrap();
        SeqAnnotation::new(
            "T".to_string(),
            "chr1".to_string(),
            Interval::new(exons[0].0, exons[exons.len() - 1].1).unwrap(),
            strand,
            exons.iter().map(iv).collect(),
            cds.iter().map(|c| CDS::new(iv(c), Phase::Zero)).collect(),
        )
    }

    /// GERP score `score(pos)` (1-based) at every needed site of `annotation`.
    fn scores(annotation: &SeqAnnotation, score: impl Fn(usize) -> Option<f32>) -> SiteValues<f32> {
        let needed = NeededSites::from_annotations(std::slice::from_ref(annotation), false);
        let mut text = String::new();
        for pos in 1..2000 {
            if needed.contains("chr1", pos - 1) {
                if let Some(s) = score(pos) {
                    text.push_str(&format!("chr1\t{}\t{}\n", pos, s));
                }
            }
        }
        SiteValues::read(text.as_bytes(), &needed, "gerp", |x: &f32| x.is_finite()).unwrap()
    }

    fn end_trunc(loftee: &Loftee, a: &SeqAnnotation, pos1: usize) -> (i64, i64, Option<f64>, bool) {
        let x = loftee
            .transcript(a)
            .assess(pos1 - 1, MutationType::Nonsense)
            .unwrap();
        assert_eq!(
            x,
            loftee
                .transcript(a)
                .assess(pos1 - 1, MutationType::FrameshiftIndel)
                .unwrap()
        );
        (
            x.bp_dist.unwrap(),
            x.dist_from_last_exon.unwrap(),
            x.gerp_dist,
            x.filters.contains(Filters::END_TRUNC),
        )
    }

    #[test]
    fn end_trunc_plus_strand() {
        // exons 101-200, 301-400, 501-600 (1-based); coding 151-560 (stop codon 558-560)
        let a = annotation(
            Strand::Plus,
            &[(100, 200), (300, 400), (500, 600)],
            &[(150, 200), (300, 400), (500, 560)],
        );
        let loftee = Loftee::new(Some(scores(&a, |_| Some(1.0))), 0.0, false);
        // last coding exon: 560 - 550 walked, coding length of the last exon 560 - 501
        assert_eq!(end_trunc(&loftee, &a, 550), (10, 10 - 59, Some(11.0), false));
        // penultimate exon: 400 - 390 + (560 - 501) walked
        assert_eq!(end_trunc(&loftee, &a, 390), (69, 10, Some(11.0 + 60.0), false));
        assert_eq!(end_trunc(&loftee, &a, 340).1, 60);
        // a cutoff above the sums: END_TRUNC within the 50 bp rule only
        let strict = Loftee::new(Some(scores(&a, |_| Some(1.0))), 100.0, false);
        assert!(end_trunc(&strict, &a, 550).3);
        assert!(end_trunc(&strict, &a, 390).3);
        assert!(!end_trunc(&strict, &a, 340).3);
        // poorly conserved last exon: END_TRUNC with cutoff 0
        let poor = Loftee::new(
            Some(scores(&a, |p| Some(if p > 500 { -2.0 } else { 1.0 }))),
            0.0,
            false,
        );
        let (_, _, gerp, filtered) = end_trunc(&poor, &a, 550);
        assert_eq!((gerp, filtered), (Some(-22.0), true));
        let x = poor.transcript(&a).assess(549, MutationType::Nonsense).unwrap();
        assert_eq!(
            (x.label().as_str(), x.is_high_confidence()),
            ("LC:END_TRUNC", false)
        );
        assert_eq!(
            x.info(),
            "GERP_DIST:-22.000,BP_DIST:10,DIST_FROM_LAST_EXON:-49,50_BP_RULE:FAIL"
        );
        assert_eq!(
            poor.counted_type(MutationType::Nonsense, Some(&x)),
            MutationType::Nonsense
        );
        let filtering = Loftee::new(Some(scores(&a, |_| Some(-1.0))), 0.0, true);
        assert_eq!(
            filtering.counted_type(MutationType::Nonsense, Some(&x)),
            MutationType::Unknown
        );
        // without scores: no END_TRUNC, no GERP_DIST
        let none = Loftee::new(None, 0.0, true);
        assert_eq!(end_trunc(&none, &a, 550), (10, -49, None, false));
    }

    #[test]
    fn end_trunc_counts_3utr_exons() {
        // a 3' UTR exon 701-800 after the last coding exon adds 99 to the distance
        // (here d = 99 - (pos - 501))
        let a = annotation(
            Strand::Plus,
            &[(100, 200), (300, 400), (500, 600), (700, 800)],
            &[(150, 200), (300, 400), (500, 560)],
        );
        let loftee = Loftee::new(Some(scores(&a, |_| Some(1.0))), 200.0, false);
        assert_eq!(
            end_trunc(&loftee, &a, 550),
            (10 + 99, 50, Some(11.0 + 100.0), true)
        );
        assert_eq!(
            end_trunc(&loftee, &a, 549),
            (11 + 99, 51, Some(12.0 + 100.0), false)
        );
        assert_eq!(end_trunc(&loftee, &a, 390).0, 10 + 59 + 99);
    }

    #[test]
    fn end_trunc_minus_strand() {
        // exons 501-600, 301-400, 101-200 in transcript order; coding 560 down to 151
        let a = annotation(
            Strand::Minus,
            &[(100, 200), (300, 400), (500, 600)],
            &[(150, 200), (300, 400), (500, 560)],
        );
        let loftee = Loftee::new(Some(scores(&a, |_| Some(1.0))), 0.0, false);
        // last coding exon: 160 - 151 walked; its coding length 200 - 151
        assert_eq!(end_trunc(&loftee, &a, 160), (9, 9 - 49, Some(10.0), false));
        // 310 - 301 in its exon and 200 - 151 in the last coding exon
        assert_eq!(end_trunc(&loftee, &a, 310), (58, 9, Some(10.0 + 50.0), false));
        assert_eq!(end_trunc(&loftee, &a, 360).1, 59);
        // in the first exon: the second exon is added whole (400 - 301)
        assert_eq!(
            end_trunc(&loftee, &a, 520),
            (19 + 99 + 49, 19 + 99, Some(20.0 + 100.0 + 50.0), false)
        );
    }

    #[test]
    fn gerp_sums_extrapolate_missing_scores() {
        let a = annotation(Strand::Plus, &[(100, 200)], &[(100, 160)]);
        // scores at odd positions only: the mean of the scored bases times the length
        let loftee = Loftee::new(
            Some(scores(&a, |p| if p % 2 == 1 { Some(p as f32) } else { None })),
            0.0,
            false,
        );
        // bases 151..=160: scored 151, 153, ..., 159 (mean 155) -> 10 * 155
        assert_eq!(end_trunc(&loftee, &a, 151).2, Some(1550.0));
        // no scores at all: a sum of 0, which is at most the cutoff 0
        let unscored = Loftee::new(Some(scores(&a, |_| None)), 0.0, false);
        let (dist, d, gerp, filtered) = end_trunc(&unscored, &a, 151);
        assert_eq!((dist, d, gerp, filtered), (9, 9 - 59, Some(0.0), true));
    }

    #[test]
    fn splice_filters() {
        // 5' UTR intron 201-300, coding 321-560, small 3' UTR intron 601-610
        let a = annotation(
            Strand::Plus,
            &[(100, 200), (300, 400), (500, 600), (610, 700)],
            &[(320, 400), (500, 560)],
        );
        let loftee = Loftee::new(None, 0.0, false);
        let t = loftee.transcript(&a);
        let label = |pos0: usize| t.assess(pos0, MutationType::SpliceSite).unwrap().label();
        assert_eq!(label(200), "LC:5UTR_SPLICE");
        assert_eq!(label(299), "LC:5UTR_SPLICE");
        assert_eq!(label(400), "HC");
        assert_eq!(label(499), "HC");
        assert_eq!(label(600), "LC:SMALL_INTRON,3UTR_SPLICE");
        assert_eq!(
            t.assess(400, MutationType::SpliceSite).unwrap().info(),
            "INTRON_SIZE:100"
        );
        assert_eq!(t.assess(400, MutationType::Missense), None);
        assert_eq!(t.assess(400, MutationType::InFrameIndel), None);
        // minus strand: the intron after the stop codon is downstream
        let m = annotation(
            Strand::Minus,
            &[(100, 200), (300, 400), (500, 600), (610, 700)],
            &[(320, 400), (500, 560)],
        );
        let t = loftee.transcript(&m);
        let label = |pos0: usize| t.assess(pos0, MutationType::SpliceSite).unwrap().label();
        assert_eq!(label(200), "LC:3UTR_SPLICE");
        assert_eq!(label(400), "HC");
        assert_eq!(label(601), "LC:SMALL_INTRON,5UTR_SPLICE");
    }
}
