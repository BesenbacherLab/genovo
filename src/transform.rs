use std::convert::TryInto;
use std::io::{BufReader, BufWriter, Write};

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use mutexpect::{Interval, Phase, SeqAnnotation, Strand, CDS};

use crate::io::{get_reader, get_writer};

#[derive(Debug, Serialize, Deserialize)]
struct GFF3Record {
    seq_id: String,
    source: String,
    seq_type: String,
    start: usize,
    end: usize,
    _score: String, // might be a "." or a floating point number. But we don't care about it anyway
    strand: char,   // either + or -
    phase: char,    // 0, 1 or 2
    attributes: Option<String>,
}

pub fn transform_gff3_annotations(
    annotations_file: &str,
    filter_for_id: Option<&str>,
    required_tags: Option<Vec<&str>>,
) -> Result<Vec<SeqAnnotation>> {
    let mut result = Vec::new();

    // dummy initialization values
    let mut current_entity_name = String::new();
    let mut current_chromosome = String::new();
    let mut current_range = Interval::new(0, 1).expect("hardcoded");
    let mut current_strand = Strand::Plus;
    let mut current_exons = Vec::new();
    let mut current_cdss = Vec::new();

    let reader = get_reader(annotations_file)
        .with_context(|| format!("failed to open file {} for reading", annotations_file))?;
    let buf_reader = BufReader::new(reader);
    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(buf_reader);
    for row_result in csv_reader.deserialize() {
        let row: GFF3Record = row_result?;
        let attributes = row.attributes.context("Missing attributes in GFF3 file")?;
        if row.seq_type.as_str() == "gene" {
            continue
        }
        let transcript_type = get_attribute(&attributes, "transcript_type").context("missing transcript_type attribute")?;
        let gene_type = get_attribute(&attributes, "gene_type").context("missing gene_type attribute")?;
        //println!("attributes, {}", attributes);
        if transcript_type != "protein_coding" || gene_type != "protein_coding" {
            continue;
        }
        if let Some(ref req_tags) = required_tags {
            if let Some(tag) = get_attribute(&attributes, "tag") {
                let mut skip = false;
                let tags: Vec<&str> = tag.split(',').collect();
                for req_tag in req_tags {
                    if !tags.contains(&req_tag){
                        //println!("skip {:#?} {:#?} ", tags, req_tag);
                        skip = true;
                        continue
                    }
                }
                if skip {
                    continue
                }
            } else {
                continue
            }
        };
        //println!("attributes, {}", attributes);
        match row.seq_type.as_str() {
            "transcript" => {
                if current_entity_name != "" {
                    current_cdss.sort();
                    current_exons.sort();
                    
                    // if we have a previous transcript (or the transcript that we filter_for_id
                    let anno = SeqAnnotation::new(
                        current_entity_name.clone(),
                        current_chromosome,
                        current_range,
                        current_strand,
                        current_exons.clone(),
                        current_cdss.clone(),
                    );
                    if let Some(id) = filter_for_id {
                        if id == current_entity_name {
                            result.push(anno);
                        }
                    } else {
                        // Ignore genes with no CDS:
                        if current_cdss.len() != 0 {
                            result.push(anno);
                        }
                    }
                }
                current_entity_name = get_attribute(&attributes, "ID")
                    .context("missing ID attribute")?
                    .to_string();
                current_chromosome = row.seq_id;
                current_range = Interval::new(row.start - 1, row.end)?; // from 1-based to 0-based. End-exclusive
                current_strand = row.strand.try_into()?;
                current_exons.clear();
                current_cdss.clear();
            }
            "exon" => {
                let id = get_attribute(&attributes, "ID")
                    .context("missing ID attribute in GFF3 file")?;
                let parent = get_attribute(&attributes, "Parent")
                    .context("missing Parent attribute in GFF3 file")?;
                if parent != current_entity_name {
                    return Err(anyhow::anyhow!(
                        "The gff3 file is not an ordered tree structure: Exon {} has parent {}",
                        id,
                        parent
                    ));
                }
                current_exons.push(Interval::new(row.start - 1, row.end)?);
            }
            "CDS" => {
                let id = get_attribute(&attributes, "ID").context("missing ID attribute")?;
                let parent =
                    get_attribute(&attributes, "Parent").context("missing Parent attribute")?;
                if parent != current_entity_name {
                    return Err(anyhow::anyhow!(
                        "The gff3 file is not an ordered tree structure: CDS {} has parent {}",
                        id,
                        parent
                    ));
                }
                let phase: Phase = row
                    .phase
                    .try_into()
                    .context("CDS region without a proper phase")?;
                current_cdss.push(CDS::new(Interval::new(row.start - 1, row.end)?, phase));
            }
            _ => {}
        }
    }
    // finish off the last entry
    current_cdss.sort();
    current_exons.sort();
    if current_entity_name != "" {
        //if we have a previous transcript
        if let Some(id) = filter_for_id {
            if id != current_entity_name {
                return Ok(result);
            }
        }
        let anno = SeqAnnotation::new(
            current_entity_name,
            current_chromosome,
            current_range,
            current_strand,
            current_exons,
            current_cdss,
        );
        result.push(anno);
    }

    Ok(result)
}

fn get_attribute<'a>(attr_str: &'a str, attribute_name: &str) -> Option<&'a str> {
    for attribute in attr_str.split(';') {
        if attribute.starts_with(attribute_name) {
            let key_value: Vec<&str> = attribute.split('=').collect();
            if key_value.len() != 2 || key_value[0] != attribute_name {
                // make sure the attribute name is not just a prefix
                continue;
            }
            return Some(key_value[1]);
        }
    }
    None
}

// this is not the best place to put it semantically, but the read() function is in the other crate
// and this uses some utility functions from *this* crate.
pub fn write_sequence_annotations_to_file(
    out_path: &str,
    sequence_annotations: &[SeqAnnotation],
) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut buf_writer = BufWriter::new(writer);
    for anno in sequence_annotations {
        let mut exons_str = String::with_capacity(anno.exons.len() * 10); // 10 is a guestimate
        for exon in &anno.exons {
            if !exons_str.is_empty() {
                exons_str.push(';')
            }
            exons_str.push_str(&exon.to_string())
        }
        let mut cds_str = String::with_capacity(anno.coding_sequences.len() * 10);
        let mut cds_phase_str = String::with_capacity(anno.coding_sequences.len() * 2); // that should be accurate
        for cds in &anno.coding_sequences {
            if !cds_str.is_empty() {
                cds_str.push(';');
                cds_phase_str.push(';');
            }
            cds_str.push_str(&cds.range.to_string());
            cds_phase_str.push(cds.phase.into());
        }
        buf_writer.write_all(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                anno.name, anno.chr, anno.strand, anno.range, exons_str, cds_str, cds_phase_str
            )
            .as_bytes(),
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gff3_io() {
        let file = "/tmp/unit_test.transform_annotations.gff3";
        let anno_str = "# some coment
# another comment
chr1	test	gene	1	100	.	+	.	attrs
chr1	test	gene	1	100	.	+	.	attrs
chr1	test	gene	1	100	.	+	.	attrs
chr1	test	transcript	10	90	.	+	.	foo=bar;ID=transcript1;baz=quux;gene_type=protein_coding;transcript_type=protein_coding
chr1	test	exon	20	30	.	+	.	ID=ex1;Parent=transcript1;gene_type=protein_coding;transcript_type=protein_coding
chr1	test	exon	35	40	.	+	.	ID=ex2;Parent=transcript1;gene_type=protein_coding;transcript_type=protein_coding
chr1	test	CDS	20	25	.	+	2	Parent=transcript1;ID=cds1;bla=bla;gene_type=protein_coding;transcript_type=protein_coding
chr1	test	CDS	38	40	.	+	1	bla=bla;Parent=transcript1;ID=cds2;gene_type=protein_coding;transcript_type=protein_coding
chr2	test	transcript	10	90	.	+	.	foo=bar;ID=transcript2;baz=quux;gene_type=protein_coding;transcript_type=protein_coding
chr2	test	exon	20	30	.	+	.	ID=ex3;Parent=transcript2;gene_type=protein_coding;transcript_type=protein_coding
chr2	test	exon	35	40	.	+	.	ID=ex4;Parent=transcript2;gene_type=protein_coding;transcript_type=protein_coding
chr2	test	CDS	20	25	.	+	2	Parent=transcript2;ID=cds2;bla=bla;gene_type=protein_coding;transcript_type=protein_coding
chr2	test	CDS	38	40	.	+	1	bla=bla;Parent=transcript2;ID=cds3;gene_type=protein_coding;transcript_type=protein_coding
";
        let mut fd = std::fs::File::create(file).unwrap();
        fd.write_all(anno_str.as_bytes()).unwrap();
        drop(fd);

        let annos = transform_gff3_annotations(file, None, None).unwrap();
        assert_eq!(annos.len(), 2);
        let a = &annos[0];
        assert_eq!(a.name, "transcript1");
        assert_eq!(a.chr, "chr1");
        assert_eq!(a.range, Interval::new(9, 90).unwrap());
        assert_eq!(a.exons.len(), 2);
        assert_eq!(a.exons[0], Interval::new(19, 30).unwrap());
        assert_eq!(a.exons[1], Interval::new(34, 40).unwrap());
        assert_eq!(a.coding_sequences.len(), 2);
        assert_eq!(a.coding_sequences[0].range, Interval::new(19, 25).unwrap());
        assert_eq!(a.coding_sequences[0].phase, Phase::Two);
        assert_eq!(a.coding_sequences[1].range, Interval::new(37, 40).unwrap());
        assert_eq!(a.coding_sequences[1].phase, Phase::One);

        let a = &annos[1];
        assert_eq!(a.name, "transcript2");
        assert_eq!(a.chr, "chr2");
        assert_eq!(a.range, Interval::new(9, 90).unwrap());
        assert_eq!(a.exons.len(), 2);
        assert_eq!(a.exons[0], Interval::new(19, 30).unwrap());
        assert_eq!(a.exons[1], Interval::new(34, 40).unwrap());
        assert_eq!(a.coding_sequences.len(), 2);
        assert_eq!(a.coding_sequences[0].range, Interval::new(19, 25).unwrap());
        assert_eq!(a.coding_sequences[0].phase, Phase::Two);
        assert_eq!(a.coding_sequences[1].range, Interval::new(37, 40).unwrap());
        assert_eq!(a.coding_sequences[1].phase, Phase::One);
    }

    #[test]
    fn test_region_file_io() {
        let file = "/tmp/unit_test.transform_annotations.regions";

        let annos = vec![
            SeqAnnotation::new(
                "transcript1".to_string(),
                "chr1".to_string(),
                Interval::new(1, 100).unwrap(),
                Strand::Plus,
                vec![
                    // exons
                    Interval::new(10, 20).unwrap(),
                    Interval::new(30, 40).unwrap(),
                ],
                vec![
                    // CDS
                    CDS::new(Interval::new(15, 20).unwrap(), Phase::Two),
                    CDS::new(Interval::new(32, 40).unwrap(), Phase::One),
                ],
            ),
            SeqAnnotation::new(
                "transcript2".to_string(),
                "chr2".to_string(),
                Interval::new(1, 100).unwrap(),
                Strand::Plus,
                vec![
                    // exons
                    Interval::new(10, 20).unwrap(),
                    Interval::new(30, 40).unwrap(),
                ],
                vec![
                    // CDS
                    CDS::new(Interval::new(15, 20).unwrap(), Phase::Two),
                    CDS::new(Interval::new(32, 40).unwrap(), Phase::One),
                ],
            ),
        ];

        write_sequence_annotations_to_file(file, &annos).unwrap();
        let annos2 = mutexpect::read_sequence_annotations_from_file(file, None).unwrap();
        assert_eq!(annos, annos2);
    }
}
