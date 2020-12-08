use std::cmp::Ordering::Equal;
use std::collections::HashMap;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::counts::{ExpectedMutationCounts, ObservedMutationCounts};
use crate::io::get_writer;
use crate::observed::AnnotatedPointMutation;
use crate::sample::SampledMutations;
use crate::{Float, MutationType};

pub fn compare_mutations(
    classified_observed_mutations: &[AnnotatedPointMutation],
    expected_mutations: &HashMap<String, ExpectedMutationCounts>,
    sampled_mutations: &SampledMutations,
    filter_for_id: Option<&str>,
) -> Result<Vec<ComparedMutations>> {
    let mut result = vec![];
    let observed_mutations =
        tally_up_observed_mutations(classified_observed_mutations, filter_for_id);

    let no_observations = ObservedMutationCounts::default(); // no observed mutations
    for (region, region_expected) in expected_mutations {
        let region_observed = observed_mutations.get(region).unwrap_or(&no_observations);
        let region_sampled = sampled_mutations.get(region).with_context(|| {
            format!("Failed to look up sampled mutations for region {}", &region)
        })?;
        for mutation_type in MutationType::iter() {
            if mutation_type == MutationType::Unknown {
                continue; // we don't sample these
            }
            let expected = region_expected.get(mutation_type);
            let observed = region_observed.get(mutation_type);
            let sampled = match region_sampled.get(&mutation_type) {
                Some(v) => v,
                None => {
                    if !expected.eq(&0.0) {
                        eprintln!("[WARNING] mutation_type={} has no samplings but an expectation value of {} in region {}", mutation_type, expected, region);
                    }
                    continue; // will not add to result
                }
            };
            let p_value = sampled.p_values().n_hits_or_more(observed);
            let comparison =
                ComparedMutations::new(region.clone(), mutation_type, observed, expected, p_value);
            result.push(comparison);
        }
    }
    result.sort_unstable_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap_or(Equal));
    Ok(result)
}

pub fn tally_up_observed_mutations(
    mutations: &[AnnotatedPointMutation],
    filter_for_id: Option<&str>,
) -> HashMap<String, ObservedMutationCounts> {
    let mut result = HashMap::new();
    for mutation in mutations {
        if let Some(id) = filter_for_id {
            if mutation.region_name != id {
                continue;
            }
        }
        result
            .entry(mutation.region_name.clone())
            .or_insert_with(ObservedMutationCounts::default)
            .add(mutation.mutation_type, 1);
    }
    result
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ComparedMutations {
    region: String,
    mutation_type: &'static str, // I want a string representation in the output file
    observed: usize,
    expected: Float,
    p_value: Float,
}

impl ComparedMutations {
    pub fn new(
        region: String,
        mutation_type: MutationType,
        observed: usize,
        expected: Float,
        p_value: Float,
    ) -> Self {
        Self {
            region,
            observed,
            expected,
            p_value,
            mutation_type: mutation_type.as_str(),
        }
    }
}

pub fn write_to_file(out_path: &str, compared_mutations: &[ComparedMutations]) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    for comparison in compared_mutations {
        csv_writer.serialize(comparison)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {}
