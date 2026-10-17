use std::collections::hash_map::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};

use anyhow::{Context, Result};

use mutexpect::{
    for_each_possible_mutation, possible_mutations_by_position, Alteration, MutationEvent,
    SeqAnnotation,
};
use pattern_partition_prediction::{PaPaPred, PaPaPredIndel};
use twobit::TwoBitFile;

use crate::error::ParseError;
use crate::io::{get_reader, get_writer};
use crate::scaling::{SiteScaling, VariantClass};
use crate::{Float, MutationType};

type PossibleMutations = HashMap<String, Vec<MutationEvent>>;

pub fn enumerate_possible_mutations(
    annotations: &[SeqAnnotation],
    ref_genome: &TwoBitFile,
    mutation_rates: &PaPaPred,
    indel_mutation_rates: &Option<PaPaPredIndel>,
    scaling: &SiteScaling,
    drop_nan: bool,
    filter_for_id: Option<&str>,
    include_intronic : bool,
    include_unknown : bool,
    filter_plof : bool,
) -> Result<PossibleMutations> {
    let mut result = HashMap::new();
    let mut masked = 0usize;
    let radius = mutation_rates.radius();
    for annotation in annotations {
        if let Some(id) = filter_for_id {
            if annotation.name != id {
                continue;
            }
        }
        let seq = match annotation_sequence(ref_genome, annotation, radius) {
            Ok(seq) => seq,
            Err(e) => {
                eprintln!(
                    "[WARNING] Skipping faulty annotation {}: {}",
                    annotation.name, e
                );
                continue;
            }
        };
        match possible_mutations_by_position(&seq, &annotation, mutation_rates, indel_mutation_rates, drop_nan, include_intronic, include_unknown, filter_plof) {
            Ok(events) => {
                let mut mutations = Vec::with_capacity(events.len());
                for (position, mut event) in events {
                    let class = VariantClass::from(event.mutation_type);
                    match scaling.probability(event.probability as f64, class, &annotation.chr, position) {
                        Some(p) => {
                            event.probability = p as f32;
                            mutations.push(event);
                        }
                        None => masked += 1,
                    }
                }
                result.insert(annotation.name.clone(), mutations);
            }
            Err(e) => {
                eprintln!(
                    "[WARNING] Skipping faulty annotation {}: {}",
                    annotation.name, e
                );
            }
        }
    }

    if masked > 0 {
        eprintln!(
            "[INFO] {} possible mutations at masked sites (allele number 0 or below the minimum) left out",
            masked
        );
    }
    Ok(result)
}

/// The genomic sequence of an annotated region plus `radius` flanking bases on
/// each side.
fn annotation_sequence(
    ref_genome: &TwoBitFile,
    annotation: &SeqAnnotation,
    radius: usize,
) -> Result<String> {
    let start = annotation.range.start.checked_sub(radius).with_context(|| {
        format!(
            "region starts less than {} bases from the start of {}",
            radius, annotation.chr
        )
    })?;
    let stop = annotation.range.stop + radius + 1;
    Ok(ref_genome.sequence(&annotation.chr, start, stop)?)
}

/// The mutation rate models and per-site scaling used to add a probability
/// column to [`write_possible_mutations`].
pub struct SiteProbabilities<'a> {
    pub point: &'a PaPaPred,
    pub indel: &'a Option<PaPaPredIndel>,
    pub scaling: &'a SiteScaling,
}

/// Write every possible mutation of the annotated regions as a tab-separated
/// line: chromosome, 1-based position, reference base, alternative, transcript
/// ID, mutation type and the sequence context of `2 * radius + 1` bases.
///
/// With `probabilities`, three more columns follow: the probability of the
/// mutation as in `enumerate` (after `--scaling-factor` or `--calibration`
/// and the allele numbers), the allele number of its site (`NA` without
/// allele numbers) and the per-generation rate of the model before scaling;
/// and the mutations that `enumerate` leaves out are not written: those at masked sites, with an undefined (NaN) rate, and
/// indels without an indel model. Summing the probabilities per transcript and
/// mutation type gives the expected numbers of `expect`.
///
/// Point mutations: the position and reference base of the mutated site, and
/// the new base. Indels are written as in VCF and in the observed-mutations
/// input: the position and base of the anchor (the base before the
/// insertion/deletion point), and as alternative the anchor followed by `NNN`
/// (in-frame) or `N` (frameshift). The context is centred on the base after
/// the insertion/deletion point.
pub fn write_possible_mutations<W: Write>(
    mut out: W,
    annotations: &[SeqAnnotation],
    ref_genome: &TwoBitFile,
    radius: usize,
    filter_for_id: Option<&str>,
    include_intronic: bool,
    include_unknown: bool,
    filter_plof: bool,
    probabilities: Option<&SiteProbabilities>,
) -> Result<()> {
    let mut masked = 0usize;
    for annotation in annotations {
        if let Some(id) = filter_for_id {
            if annotation.name != id {
                continue;
            }
        }
        let seq = match annotation_sequence(ref_genome, annotation, radius) {
            Ok(seq) => seq,
            Err(e) => {
                eprintln!(
                    "[WARNING] Skipping faulty annotation {}: {}",
                    annotation.name, e
                );
                continue;
            }
        };
        // Lines of this annotation; written only if all its sequence contexts
        // have rates. Like enumerate (and calibrate), an annotation with a
        // context that cannot be looked up (e.g. an N in the reference) is
        // skipped as a whole, so the probabilities still sum to `expect`.
        let mut lines: Vec<u8> = Vec::new();
        let mut annotation_masked = 0usize;
        let result = for_each_possible_mutation(
            &seq,
            annotation,
            radius,
            include_intronic,
            include_unknown,
            filter_plof,
            |site| -> Result<()> {
                let out = &mut lines;
                let point_rates = match probabilities {
                    Some(p) => Some(p.point.rates(site.seq_window).with_context(|| {
                        format!("bad sequence {} in {}", site.seq_window, annotation.name)
                    })?),
                    None => None,
                };
                let mut indel_rates = None;
                for m in site.mutations {
                    let probability = match (probabilities, &point_rates) {
                        (Some(p), Some(rates)) => {
                            let rate = match m.alteration {
                                Alteration::Point(base) => rates.by_name(base.name()),
                                Alteration::InFrameIndel | Alteration::FrameshiftIndel => {
                                    let model = match p.indel {
                                        Some(model) => model,
                                        None => continue, // as enumerate
                                    };
                                    if indel_rates.is_none() {
                                        let w = site.seq_window;
                                        indel_rates = Some(model.rates(&w[0..w.len() - 1]).with_context(|| {
                                            format!("bad sequence {} in {}", w, annotation.name)
                                        })?);
                                    }
                                    let r = indel_rates.as_ref().expect("set above");
                                    if m.alteration == Alteration::InFrameIndel {
                                        r.inframe
                                    } else {
                                        r.outframe
                                    }
                                }
                            };
                            if rate.is_nan() {
                                continue; // as enumerate (drop_nan)
                            }
                            let class = VariantClass::from(m.mutation_type);
                            match p.scaling.probability(rate as f64, class, &annotation.chr, site.genomic_position) {
                                Some(probability) => {
                                    let an = p.scaling.allele_number(&annotation.chr, site.genomic_position);
                                    Some((probability as f32, an, rate))
                                }
                                None => {
                                    annotation_masked += 1;
                                    continue;
                                }
                            }
                        }
                        _ => None,
                    };
                    // 1-based position and reference base of the written site
                    let (position, ref_base, alt) = match m.alteration {
                        Alteration::Point(base) => (
                            site.genomic_position + 1,
                            site.ref_base.name(),
                            base.name().to_string(),
                        ),
                        Alteration::InFrameIndel | Alteration::FrameshiftIndel => {
                            // the indel lies between genomic_position - 1 and
                            // genomic_position; the anchor is the base before
                            let anchor = site.seq_window.as_bytes()[radius - 1] as char;
                            let n = if m.alteration == Alteration::InFrameIndel {
                                "NNN"
                            } else {
                                "N"
                            };
                            (site.genomic_position, anchor, format!("{}{}", anchor, n))
                        }
                    };
                    write!(
                        out,
                        "{}\t{}\t{}\t{}\t{}\t{:?}\t{}",
                        annotation.chr,
                        position,
                        ref_base,
                        alt,
                        annotation.name,
                        m.mutation_type,
                        site.seq_window
                    )?;
                    match probability {
                        Some((p, Some(an), rate)) => writeln!(out, "\t{:e}\t{}\t{:e}", p, an, rate)?,
                        Some((p, None, rate)) => writeln!(out, "\t{:e}\tNA\t{:e}", p, rate)?,
                        None => writeln!(out)?,
                    }
                }
                Ok(())
            },
        );
        match result {
            Ok(()) => {
                out.write_all(&lines)?;
                masked += annotation_masked;
            }
            Err(e) => eprintln!("[WARNING] Skipping annotation {}: {:#}", annotation.name, e),
        }
    }
    out.flush()?;
    if masked > 0 {
        eprintln!(
            "[INFO] {} possible mutations at masked sites (allele number 0 or below the minimum) left out",
            masked
        );
    }
    Ok(())
}

pub fn write_to_file(out_path: &str, possible_mutations: &PossibleMutations) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut buf_writer = BufWriter::new(writer);
    let mut names: Vec<&String> = possible_mutations.keys().collect();
    names.sort();
    for name in names {
        let mutations = &possible_mutations[name];
        buf_writer.write_all(format!("#{}\n", name).as_bytes())?;
        for mutation in mutations {
            buf_writer.write_all(
                format!(
                    "{}:{}\n",
                    mutation.mutation_type as usize, mutation.probability
                )
                .as_bytes(),
            )?;
        }
    }

    Ok(())
}

pub fn read_from_file(in_path: &str) -> Result<PossibleMutations> {
    let mut result: PossibleMutations = HashMap::new();
    let reader = get_reader(in_path)
        .with_context(|| format!("failed to open file {} for reading", in_path))?;
    let bufreader = BufReader::new(reader);
    let mut current_gene: Option<String> = None;
    for (line_no, line) in bufreader.lines().enumerate() {
        if let Ok(line) = line {
            if line.starts_with('#') {
                current_gene = Some(line[1..].to_string());
                result.insert(line[1..].to_string(), vec![]);
                continue;
            }
            if current_gene.is_none() {
                return Err(ParseError::new(format!(
                    "Expected #name line on line {} in file {}",
                    line_no + 1,
                    in_path
                ))
                .into());
            }
            let tokens: Vec<&str> = line.split(':').collect();
            let mut_type: MutationType = tokens[0].parse::<u8>()?.into();
            let probability: Float = tokens[1].parse()?;
            if let Some(gene) = &current_gene {
                result
                    .get_mut(gene)
                    .expect("init")
                    .push(MutationEvent::new(mut_type, probability));
            //result.entry(gene.clone()).or_insert_with(||Vec::new()).push(
            } else {
                panic!("gene must be Some because of previous check");
            }
        } else {
            Err(line.unwrap_err()).with_context(|| format!("IO error for file {}", in_path))?
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::convert::TryInto;

    #[test]
    fn test_possible_mutations_io() {
        use mutexpect::MutationEvent;

        //external constructor
        fn mevent(mut_type: &str, probability: Float) -> MutationEvent {
            MutationEvent::new(mut_type.try_into().unwrap(), probability)
        }

        //this is not portable and will not work on Windows
        let path = "/tmp/unit_test.possible_mutations";
        let mut pm: PossibleMutations = HashMap::new();
        write_to_file(path, &pm).unwrap();
        let pm2 = read_from_file(path).unwrap();
        assert_eq!(pm, pm2);

        pm.insert("foo".to_string(), vec![]);
        write_to_file(path, &pm).unwrap();
        let pm2 = read_from_file(path).unwrap();
        assert_eq!(pm, pm2);

        let events = vec![
            mevent("synonymous", 0.1),
            mevent("missense", 0.2),
            mevent("nonsense", 0.3),
        ];
        pm.insert("bar".to_string(), events);

        write_to_file(path, &pm).unwrap();
        let pm2 = read_from_file(path).unwrap();
        assert_eq!(pm, pm2);

        pm.insert("baz".to_string(), vec![]);
        write_to_file(path, &pm).unwrap();
        let pm2 = read_from_file(path).unwrap();
        assert_eq!(pm, pm2);
    }
}
