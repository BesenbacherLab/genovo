use std::collections::hash_map::HashMap;
use std::str::FromStr;
use std::convert::TryInto;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use bigdecimal::{BigDecimal, ToPrimitive};

use mutexpect::{MutationEvent, MutationType};

use crate::counts::ExpectedMutationCounts;
use crate::io::{get_reader, get_writer};
use crate::Float;

pub fn expected_number_of_mutations(
    possible_mutations: &HashMap<String, Vec<MutationEvent>>,
    filter_for_id: Option<&str>,
) -> Result<HashMap<String, ExpectedMutationCounts>> {
    let mut result = HashMap::new();
    for (region, events) in possible_mutations {
        if let Some(id) = filter_for_id {
            if region != id {
                continue;
            }
        }

        let mut counts = ExpectedMutationCounts::default();
        let mut precise_counts: HashMap<MutationType, BigDecimal> = HashMap::new();
        let zero = BigDecimal::from_str("0.0").expect("trivial");
        for event in events {
            //counts.add(event.mutation_type, event.probability as Float); //TODO remove me

            if !precise_counts.contains_key(&event.mutation_type) {
                precise_counts.insert(event.mutation_type, zero.clone());
            }
            let value: BigDecimal = event.probability.try_into().unwrap_or(zero.clone());
            *(precise_counts.get_mut(&event.mutation_type).unwrap()) += value;

        }

        for (mutation_type, probability) in precise_counts {
            counts.add(mutation_type, probability.to_f32().unwrap());
        }
        result.insert(region.clone(), counts);

    }
    Ok(result)
}

// serialization stuff //

type Count = Float;
#[derive(Debug, Serialize, Deserialize)]
struct CSVRow {
    region: String,
    unknown: Count,
    synonymous: Count,
    missense: Count,
    nonsense: Count,
    start_codon: Count,
    stop_loss: Count,
    splice_site: Count,
    intronic: Count,
    inframe_indel: Count,
    frameshift_indel: Count,
}

pub fn write_to_file(
    out_path: &str,
    expected_mutations: &HashMap<String, ExpectedMutationCounts>,
) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    let mut regions: Vec<&String> = expected_mutations.keys().collect();
    regions.sort();
    for region in regions {
        let c = &expected_mutations[region];
        let row = CSVRow {
            region: region.clone(),
            unknown: c.unknown,
            synonymous: c.synonymous,
            missense: c.missense,
            nonsense: c.nonsense,
            start_codon: c.start_codon,
            stop_loss: c.stop_loss,
            splice_site: c.splice_site,
            intronic: c.intronic,
            inframe_indel: c.inframe_indel,
            frameshift_indel: c.frameshift_indel,
        };
        csv_writer.serialize(row)?;
    }
    Ok(())
}

pub fn read_from_file(in_path: &str) -> Result<HashMap<String, ExpectedMutationCounts>> {
    let mut result = HashMap::new();
    let reader = get_reader(in_path)
        .with_context(|| format!("failed to open file {} for reading", in_path))?;
    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    for row_result in csv_reader.deserialize() {
        let row: CSVRow = row_result?;
        let counts = ExpectedMutationCounts {
            unknown: row.unknown,
            synonymous: row.synonymous,
            missense: row.missense,
            nonsense: row.nonsense,
            start_codon: row.start_codon,
            stop_loss: row.stop_loss,
            splice_site: row.splice_site,
            intronic: row.intronic,
            inframe_indel: row.inframe_indel,
            frameshift_indel: row.frameshift_indel,
        };
        result.insert(row.region, counts);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expected_mutations_io() {
        let path = "/tmp/unit_test.expected_mutations";
        let mut em: HashMap<String, ExpectedMutationCounts> = HashMap::new();
        write_to_file(path, &em).unwrap();
        let em2 = read_from_file(path).unwrap();
        assert_eq!(em2.len(), 0);

        let counts1 = ExpectedMutationCounts {
            unknown: 1.2,
            synonymous: 2.3,
            missense: 3.4,
            nonsense: 4.5,
            start_codon: 5.6,
            stop_loss: 6.7,
            splice_site: 7.8,
            intronic: 8.9,
            inframe_indel: 9.1,
            frameshift_indel: 10.11,
        };
        em.insert("foo".to_string(), counts1);
        write_to_file(path, &em).unwrap();
        let em2 = read_from_file(path).unwrap();
        assert_eq!(em2.len(), 1);
        assert_eq!(em2["foo"], counts1);

        let counts2 = ExpectedMutationCounts {
            unknown: 2.2,
            synonymous: 3.3,
            missense: 4.4,
            nonsense: 5.5,
            start_codon: 6.6,
            stop_loss: 7.7,
            splice_site: 8.8,
            intronic: 9.9,
            inframe_indel: 10.10,
            frameshift_indel: 11.11,
        };
        em.insert("bar".to_string(), counts2);
        write_to_file(path, &em).unwrap();
        let em2 = read_from_file(path).unwrap();
        assert_eq!(em2.len(), 2);
        assert_eq!(em2["foo"], counts1);
        assert_eq!(em2["bar"], counts2);
    }
}
