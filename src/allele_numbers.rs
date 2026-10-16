//! Per-site allele numbers: the number of called chromosomes at each site of
//! a population data set (see `docs/design/allele-number-scaling.md`).
//!
//! The input is one tab-separated file, sorted by chromosome and position,
//! with 1-based positions, either as `chrom  pos  an` or as `chrom:pos  an`
//! (the gnomAD all-sites AN layout). An optional header line and lines
//! starting with `#` are skipped. Sites that are not in the file have an
//! allele number of 0 (not callable).
//!
//! The file may cover the whole genome, so it is streamed once and only the
//! positions genovo needs (see [`NeededSites`]) are kept.
//!
//! Step 1 of the implementation plan: the reader is used by `enumerate` and
//! `classify` from step 3 on.
#![allow(dead_code)]

use std::collections::HashMap;
use std::convert::TryFrom;
use std::io::{BufRead, BufReader, Write};

use anyhow::{bail, Context, Result};

use mutexpect::interval::Interval;
use mutexpect::SeqAnnotation;

use crate::io::{get_reader, get_writer};

/// Bases kept on each side of every exon: covers the splice sites (2 intronic
/// bases), the anchor base of indels just before an exon, and the extra
/// position genovo visits after the end of a region.
pub const EXON_MARGIN: usize = 2;

/// The genomic intervals (0-based, half-open, sorted and merged per
/// chromosome) where allele numbers are needed.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct NeededSites {
    /// Chromosomes in order of first appearance in the annotations.
    chromosomes: Vec<String>,
    intervals: HashMap<String, Vec<Interval>>,
}

impl NeededSites {
    /// Exons ± [`EXON_MARGIN`] of all annotations, or the whole annotated
    /// regions (± the margin) if intronic sites are needed too.
    pub fn from_annotations(annotations: &[SeqAnnotation], include_intronic: bool) -> Self {
        let mut raw: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
        let mut chromosomes = Vec::new();
        for annotation in annotations {
            let list = raw.entry(annotation.chr.clone()).or_insert_with(|| {
                chromosomes.push(annotation.chr.clone());
                Vec::new()
            });
            let mut add = |range: &Interval| {
                list.push((
                    range.start.saturating_sub(EXON_MARGIN),
                    range.stop + EXON_MARGIN,
                ))
            };
            if include_intronic {
                add(&annotation.range);
            } else {
                annotation.exons.iter().for_each(&mut add);
                // a CDS outside the listed exons should not happen, but be safe
                annotation
                    .coding_sequences
                    .iter()
                    .for_each(|c| add(&c.range));
            }
        }
        let intervals = raw
            .into_iter()
            .map(|(chr, mut list)| {
                list.sort_unstable();
                let mut merged: Vec<Interval> = Vec::with_capacity(list.len());
                for (start, stop) in list {
                    match merged.last_mut() {
                        Some(last) if start <= last.stop => last.stop = last.stop.max(stop),
                        _ => merged.push(Interval { start, stop }),
                    }
                }
                (chr, merged)
            })
            .collect();
        NeededSites {
            chromosomes,
            intervals,
        }
    }

    /// Total number of positions.
    pub fn len(&self) -> usize {
        self.intervals.values().flatten().map(|i| i.len()).sum()
    }

    pub fn contains(&self, chr: &str, pos: usize) -> bool {
        self.intervals
            .get(chr)
            .map_or(false, |list| find(list, pos).is_some())
    }

    /// Write the intervals as BED (0-based, half-open), e.g. for restricting a
    /// genome-wide allele number file with `tabix -R` or `bedtools intersect`.
    pub fn write_bed(&self, path: &str) -> Result<()> {
        let mut out = get_writer(path)
            .with_context(|| format!("failed to open file {} for writing", path))?;
        for chr in &self.chromosomes {
            for i in &self.intervals[chr] {
                writeln!(out, "{}\t{}\t{}", chr, i.start, i.stop)?;
            }
        }
        out.flush()?;
        Ok(())
    }
}

/// Index of the interval in the sorted, disjoint `list` containing `pos`.
fn find(list: &[Interval], pos: usize) -> Option<usize> {
    let i = list.partition_point(|iv| iv.stop <= pos);
    (i < list.len() && list[i].start <= pos).then(|| i)
}

/// Allele numbers at the needed sites.
#[derive(Debug, Default)]
pub struct AlleleNumbers {
    /// Per chromosome: (0-based position, allele number), sorted by position.
    sites: HashMap<String, Vec<(u32, u32)>>,
}

impl AlleleNumbers {
    /// The allele number at a 0-based position; 0 if the site is not in the
    /// file (not callable) or was not needed.
    pub fn get(&self, chr: &str, pos: usize) -> u32 {
        let list = match self.sites.get(chr) {
            Some(list) => list,
            None => return 0,
        };
        match list.binary_search_by_key(&(pos as u64), |&(p, _)| u64::from(p)) {
            Ok(i) => list[i].1,
            Err(_) => 0,
        }
    }

    /// Number of stored sites (with allele number > 0).
    pub fn len(&self) -> usize {
        self.sites.values().map(Vec::len).sum()
    }

    pub fn read_from_file(path: &str, needed: &NeededSites) -> Result<Self> {
        let reader = get_reader(path)
            .with_context(|| format!("failed to open file {} for reading", path))?;
        Self::read(BufReader::with_capacity(1 << 20, reader), needed)
            .with_context(|| format!("reading allele numbers from {}", path))
    }

    pub fn read<R: BufRead>(mut reader: R, needed: &NeededSites) -> Result<Self> {
        let mut sites: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
        let mut finished: Vec<String> = Vec::new(); // chromosomes already seen
        let mut current_chr = String::new();
        let mut current: Option<&Vec<Interval>> = None; // needed intervals of current_chr
        let mut next_interval = 0; // first interval that may still contain sites
        let mut last_pos: Option<usize> = None;
        let mut line = String::new();
        let mut line_no = 0;
        let mut first_data_line = true;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            line_no += 1;
            let text = line.trim_end_matches(&['\n', '\r'][..]);
            if text.is_empty() || text.starts_with('#') {
                continue;
            }
            let (chr, pos, an) = match parse_line(text) {
                Some(fields) => fields,
                // a header has the same number of columns, but no numbers
                None if first_data_line && matches!(text.split('\t').count(), 2 | 3) => {
                    first_data_line = false;
                    continue;
                }
                None => bail!(
                    "line {}: expected 'chrom<TAB>pos<TAB>an' or 'chrom:pos<TAB>an', found '{}'",
                    line_no,
                    text
                ),
            };
            first_data_line = false;
            if pos == 0 {
                bail!("line {}: positions are 1-based, found 0", line_no);
            }
            let pos = pos - 1;
            if chr != current_chr {
                if finished.iter().any(|c| c == chr) {
                    bail!(
                        "line {}: chromosome {} appears again after other chromosomes; \
                         the file must be sorted by chromosome and position",
                        line_no,
                        chr
                    );
                }
                if !current_chr.is_empty() {
                    finished.push(std::mem::take(&mut current_chr));
                }
                current_chr = chr.to_string();
                current = needed.intervals.get(chr);
                next_interval = 0;
                last_pos = None;
            }
            if let Some(prev) = last_pos {
                if pos <= prev {
                    bail!(
                        "line {}: position {} of {} is not after the previous position {}; \
                         the file must be sorted by chromosome and position",
                        line_no,
                        pos + 1,
                        chr,
                        prev + 1
                    );
                }
            }
            last_pos = Some(pos);
            let intervals = match current {
                Some(list) => list,
                None => continue, // chromosome not needed
            };
            // positions increase, so needed intervals ending before pos are done
            while next_interval < intervals.len() && intervals[next_interval].stop <= pos {
                next_interval += 1;
            }
            if next_interval < intervals.len() && intervals[next_interval].start <= pos && an > 0 {
                let pos32 = u32::try_from(pos).context("position does not fit in 32 bits")?;
                sites
                    .entry(current_chr.clone())
                    .or_insert_with(Vec::new)
                    .push((pos32, an));
            }
        }
        Ok(AlleleNumbers { sites })
    }
}

/// Parse `chrom<TAB>pos<TAB>an` or `chrom:pos<TAB>an`.
fn parse_line(text: &str) -> Option<(&str, usize, u32)> {
    let mut fields = text.split('\t');
    let first = fields.next()?;
    let second = fields.next()?;
    let (chr, pos, an) = match fields.next() {
        Some(third) => (first, second, third),
        None => {
            let (chr, pos) = first.rsplit_once(':')?;
            (chr, pos, second)
        }
    };
    if chr.is_empty() || fields.next().is_some() {
        return None;
    }
    Some((chr, pos.trim().parse().ok()?, an.trim().parse().ok()?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use mutexpect::{Phase, Strand, CDS};

    fn annotation(chr: &str, range: (usize, usize), exons: &[(usize, usize)]) -> SeqAnnotation {
        let iv = |(s, e): (usize, usize)| Interval::new(s, e).unwrap();
        SeqAnnotation::new(
            format!("T{}", range.0),
            chr.to_string(),
            iv(range),
            Strand::Plus,
            exons.iter().map(|&e| iv(e)).collect(),
            vec![CDS::new(iv(exons[0]), Phase::Zero)],
        )
    }

    fn needed() -> NeededSites {
        NeededSites::from_annotations(
            &[
                annotation("chr1", (100, 200), &[(100, 120), (180, 200)]),
                annotation("chr1", (110, 130), &[(110, 130)]), // overlaps the first
                annotation("chr2", (50, 60), &[(50, 60)]),
            ],
            false,
        )
    }

    #[test]
    fn needed_sites_are_exons_with_margin_merged() {
        let n = needed();
        assert_eq!(
            n.intervals["chr1"],
            vec![
                Interval {
                    start: 98,
                    stop: 132
                },
                Interval {
                    start: 178,
                    stop: 202
                }
            ]
        );
        assert_eq!(
            n.intervals["chr2"],
            vec![Interval {
                start: 48,
                stop: 62
            }]
        );
        assert_eq!(n.len(), 34 + 24 + 14);
        assert!(n.contains("chr1", 98) && n.contains("chr1", 131) && !n.contains("chr1", 132));
        assert!(!n.contains("chr1", 150) && !n.contains("chr3", 50));
        let whole = NeededSites::from_annotations(
            &[annotation("chr1", (100, 200), &[(100, 120), (180, 200)])],
            true,
        );
        assert_eq!(
            whole.intervals["chr1"],
            vec![Interval {
                start: 98,
                stop: 202
            }]
        );
    }

    #[test]
    fn bed_output() {
        let path = std::env::temp_dir().join("genovo_needed_sites_test.bed");
        let path = path.to_str().unwrap();
        needed().write_bed(path).unwrap();
        assert_eq!(
            std::fs::read_to_string(path).unwrap(),
            "chr1\t98\t132\nchr1\t178\t202\nchr2\t48\t62\n"
        );
    }

    fn read(text: &str) -> Result<AlleleNumbers> {
        AlleleNumbers::read(text.as_bytes(), &needed())
    }

    #[test]
    fn reads_both_layouts_keeps_needed_sites_only() {
        // 1-based positions: chr1:99 is 0-based 98, the first needed site.
        let generic = "chrom\tpos\tan\nchr1\t97\t5\nchr1\t99\t10\nchr1\t100\t12\n\
                       chr1\t150\t7\nchr1\t202\t9\nchr1\t203\t9\nchr2\t49\t3\nchr3\t1\t1\n";
        let locus = "locus\tAN\nchr1:97\t5\nchr1:99\t10\nchr1:100\t12\n\
                     chr1:150\t7\nchr1:202\t9\nchr1:203\t9\nchr2:49\t3\nchr3:1\t1\n";
        for text in [generic, locus] {
            let a = read(text).unwrap();
            assert_eq!(a.get("chr1", 98), 10);
            assert_eq!(a.get("chr1", 99), 12);
            assert_eq!(a.get("chr1", 201), 9); // last needed position
            assert_eq!(a.get("chr1", 96), 0); // not needed
            assert_eq!(a.get("chr1", 149), 0); // intronic, not needed
            assert_eq!(a.get("chr1", 202), 0);
            assert_eq!(a.get("chr2", 48), 3);
            assert_eq!(a.get("chr1", 100), 0); // needed but missing: not callable
            assert_eq!(a.get("chr3", 0), 0);
            assert_eq!(a.len(), 4);
        }
    }

    #[test]
    fn comments_without_header_and_zero_an() {
        let a = read("# a comment\nchr1\t99\t10\nchr1\t100\t0\n\nchr1\t101\t4\r\n").unwrap();
        assert_eq!(
            (a.get("chr1", 98), a.get("chr1", 99), a.get("chr1", 100)),
            (10, 0, 4)
        );
        assert_eq!(a.len(), 2); // allele number 0 is not stored
    }

    #[test]
    fn errors() {
        let unsorted = "chr1\t100\t1\nchr1\t99\t1\n";
        let msg = format!("{:#}", read(unsorted).unwrap_err());
        assert!(msg.contains("sorted"), "{}", msg);
        assert!(read("chr1\t99\t1\nchr1\t99\t2\n").is_err()); // duplicate
        let regrouped = "chr1\t99\t1\nchr2\t49\t1\nchr1\t120\t1\n";
        assert!(format!("{:#}", read(regrouped).unwrap_err()).contains("again"));
        assert!(read("chr1\t99\t1\nchr1\tx\t1\n").is_err()); // bad line after data
        assert!(read("chr1\t0\t1\n").is_err()); // not 1-based
        assert!(read("chr1\t99\t1\t5\n").is_err()); // too many columns
    }

    #[test]
    fn gzip_and_bgzip_files() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let dir = std::env::temp_dir();
        let gz = |members: &[&str]| -> Vec<u8> {
            // bgzip files are several gzip members concatenated
            let mut out = Vec::new();
            for m in members {
                let mut e = GzEncoder::new(Vec::new(), Compression::default());
                e.write_all(m.as_bytes()).unwrap();
                out.extend(e.finish().unwrap());
            }
            out
        };
        let path = dir.join("genovo_an_test.tsv.gz");
        let path = path.to_str().unwrap();
        std::fs::write(
            path,
            gz(&["chr1\t99\t10\n", "chr1\t100\t12\nchr2\t49\t3\n"]),
        )
        .unwrap();
        let a = AlleleNumbers::read_from_file(path, &needed()).unwrap();
        assert_eq!(
            (a.get("chr1", 98), a.get("chr1", 99), a.get("chr2", 48)),
            (10, 12, 3)
        );
    }
}
