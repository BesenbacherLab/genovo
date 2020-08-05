use std::collections::HashMap;
use std::convert::TryInto;
use std::default::Default;

use anyhow::{Context, Result};
use rand::Rng;
use serde::{Deserialize, Serialize};

use mutexpect::MutationEvent;

use crate::counts::DefaultCounter;
use crate::io::{get_reader, get_writer};
use crate::{Float, MutationType};

type Probability = Float;
pub type SampledMutations = HashMap<String, HashMap<MutationType, DefaultCounter>>;

pub fn sample_mutations(
    possible_mutations: &HashMap<String, Vec<MutationEvent>>,
    number_of_samplings: usize,
    drop_unknown_mutation_type: bool,
    filter_for_id: Option<&str>,
) -> Result<SampledMutations> {
    let mut random = rand::thread_rng();

    let mut result: SampledMutations = HashMap::new();
    for (region_name, mutation_events) in possible_mutations {
        if let Some(id) = filter_for_id {
            if region_name != id {
                continue;
            }
        }

        // sort mutations by mutation type
        let mut mutations: HashMap<MutationType, Vec<Probability>> = HashMap::new();
        for event in mutation_events {
            if drop_unknown_mutation_type && event.mutation_type == MutationType::Unknown {
                continue; // it makes little sense to sample these
            }
            mutations
                .entry(event.mutation_type)
                .or_insert_with(Vec::new)
                .push(event.probability);
        }

        let mut distributions = HashMap::new();

        // for each mutation type
        for (mutation_type, mutation_probabilities) in &mut mutations {
            // sort mutation probabilities
            mutation_probabilities.sort_unstable_by(|a, b| a.partial_cmp(b).expect("no NaNs")); // make binary search possible

            let mut counter = DefaultCounter::new();
            // do the random samplings
            for _ in 0..number_of_samplings {
                let exclusive_threshold = random.gen_range(0.0, 1.0); // floating point coin

                // because the values are sorted, we can simply determine the number of values that
                // are above the threshold and treat them as the number of failures.
                let failures = match mutation_probabilities.binary_search_by(|other| {
                    other.partial_cmp(&exclusive_threshold).expect("No NaNs")
                }) {
                    Ok(index) => {
                        let mut i = index; // if the same probability exists several times, we need to use the leftmost one
                        #[allow(clippy::float_cmp)]
                        while i > 1 && mutation_probabilities[i - 1] == exclusive_threshold {
                            i -= 1;
                        }
                        i // leftmost index of value
                    }
                    Err(index) => index, // minimum index of those entries higher than `threshold`
                };
                let successes = mutation_probabilities.len() - failures;
                counter.inc(successes);
            }
            distributions.insert(*mutation_type, counter);
        }
        result.insert(region_name.clone(), distributions);
    }
    Ok(result)
}

// serialization stuff //

#[derive(Debug, Serialize, Deserialize, Default)]
struct CSVRow {
    region: String,
    unknown: String,
    synonymous: String,
    missense: String,
    nonsense: String,
    start_codon: String,
    stop_loss: String,
    splice_site: String,
}

pub fn write_to_file(out_path: &str, sampled_mutations: &SampledMutations) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    for (region, distributions) in sampled_mutations {
        let mut row = CSVRow::default();
        row.region = region.clone();
        for (mut_type, distribution) in distributions {
            let value = distribution.to_string();
            match mut_type {
                MutationType::Unknown => row.unknown = value,
                MutationType::Synonymous => row.synonymous = value,
                MutationType::Missense => row.missense = value,
                MutationType::Nonsense => row.nonsense = value,
                MutationType::StartCodon => row.start_codon = value,
                MutationType::StopLoss => row.stop_loss = value,
                MutationType::SpliceSite => row.splice_site = value,
            }
        }
        csv_writer.serialize(row)?;
    }
    Ok(())
}

pub fn read_from_file(in_path: &str) -> Result<SampledMutations> {
    let mut result = SampledMutations::new();
    let reader = get_reader(in_path)
        .with_context(|| format!("failed to open file {} for reading", in_path))?;

    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    for row_result in csv_reader.deserialize() {
        let row: CSVRow = row_result?;
        let mut samplings: HashMap<MutationType, DefaultCounter> = HashMap::new();
        samplings.insert(MutationType::Unknown, row.unknown.as_str().try_into()?);
        samplings.insert(
            MutationType::Synonymous,
            row.synonymous.as_str().try_into()?,
        );
        samplings.insert(MutationType::Missense, row.missense.as_str().try_into()?);
        samplings.insert(MutationType::Nonsense, row.nonsense.as_str().try_into()?);
        samplings.insert(
            MutationType::StartCodon,
            row.start_codon.as_str().try_into()?,
        );
        samplings.insert(MutationType::StopLoss, row.stop_loss.as_str().try_into()?);
        samplings.insert(
            MutationType::SpliceSite,
            row.splice_site.as_str().try_into()?,
        );
        result.insert(row.region, samplings);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sampled_mutations_io() {
        let path = "/tmp/unit_test.sampled_mutations";
        let mut sm = SampledMutations::new();
        write_to_file(path, &sm).unwrap();
        let sm2 = read_from_file(path).unwrap();
        assert_eq!(sm2.len(), 0);

        let mut distributions1 = HashMap::new();
        distributions1.insert(MutationType::Unknown, "0|1|2|3|4".try_into().unwrap());
        distributions1.insert(MutationType::Synonymous, "10|20|30|40".try_into().unwrap());
        distributions1.insert(MutationType::Missense, "11|22|33|44".try_into().unwrap());
        distributions1.insert(MutationType::Nonsense, "2|3".try_into().unwrap());
        distributions1.insert(MutationType::StartCodon, "4".try_into().unwrap());
        distributions1.insert(MutationType::StopLoss, "5|6".try_into().unwrap());
        distributions1.insert(MutationType::SpliceSite, "0|0|5|6|7".try_into().unwrap());
        sm.insert("foo".to_string(), distributions1.clone());
        write_to_file(path, &sm).unwrap();
        let sm2 = read_from_file(path).unwrap();
        assert_eq!(sm2.len(), 1);
        assert_eq!(sm2["foo"], distributions1);

        let mut distributions2 = HashMap::new();
        distributions2.insert(MutationType::Unknown, "1|1|2|3|4".try_into().unwrap());
        distributions2.insert(MutationType::Synonymous, "11|20|30|40".try_into().unwrap());
        distributions2.insert(MutationType::Missense, "12|22|33|44".try_into().unwrap());
        distributions2.insert(MutationType::Nonsense, "6|7".try_into().unwrap());
        distributions2.insert(MutationType::StartCodon, "5".try_into().unwrap());
        distributions2.insert(MutationType::StopLoss, "6|6".try_into().unwrap());
        distributions2.insert(MutationType::SpliceSite, "0|0|6|6|7".try_into().unwrap());
        sm.insert("bar".to_string(), distributions2.clone());
        write_to_file(path, &sm).unwrap();
        let sm2 = read_from_file(path).unwrap();
        assert_eq!(sm2.len(), 2);
        assert_eq!(sm2["foo"], distributions1);
        assert_eq!(sm2["bar"], distributions2);
    }
}
