use anyhow::Result;
use clap::{App, Arg};

use pattern_partition_prediction::{PaPaPred, PaPaPredIndel};
use twobit::TwoBitFile;

mod allele_numbers;
mod calibrate;
mod compare;
mod counts;
mod enumerate;
mod error;
mod expect;
mod io;
mod observed;
mod sample;
mod scaling;
mod seed;
mod transform;

const VERSION: &'static str = env!("CARGO_PKG_VERSION");

type Float = f32;
type MutationType = mutexpect::MutationType;

use crate::allele_numbers::{AlleleNumbers, NeededSites};
use crate::scaling::{Calibration, SiteScaling};
use crate::compare::compare_mutations;
use crate::enumerate::enumerate_possible_mutations;
use crate::enumerate::write_possible_mutations;
use crate::error::MissingCommandLineArgumentError;
use crate::expect::expected_number_of_mutations;
use crate::observed::classify_mutations;
use crate::observed::read_mutations_from_file as read_observed_mutations_from_file;
use crate::sample::sample_mutations;

fn require_initialization<'a, T>(
    value: &'a Option<T>,
    cli_argument: &'static str,
) -> std::result::Result<&'a T, MissingCommandLineArgumentError> {
    match value {
        Some(v) => Ok(v),
        None => Err(MissingCommandLineArgumentError::new(cli_argument)),
    }
}

/// `--action calibrate`: fit S and h on the possible synonymous mutations.
fn run_calibrate(
    matches: &clap::ArgMatches,
    regions: &Option<Vec<mutexpect::SeqAnnotation>>,
    genome: &Option<TwoBitFile>,
    papa: &Option<PaPaPred>,
    id: Option<&str>,
) -> Result<()> {
    let regions = require_initialization(regions, "--genomic-regions or --gff3")?;
    let genome = require_initialization(genome, "--genome")?;
    let papa = require_initialization(papa, "--point-mutation-probabilities")?;
    let required = |name: &'static str| {
        matches
            .value_of(name)
            .ok_or_else(|| anyhow::anyhow!("--action calibrate needs --{}", name))
    };
    let an_file = required("allele-numbers")?;
    let observed_file = required("observed-mutations")?;
    let out = required("calibration")?;
    let n_bins: usize = matches
        .value_of("calibration-bins")
        .expect("default value")
        .parse()
        .map_err(|_| anyhow::anyhow!("--calibration-bins must be a positive integer"))?;
    let n_min: u32 = match matches.value_of("min-an") {
        Some(v) => v
            .parse()
            .map_err(|_| anyhow::anyhow!("--min-an must be a non-negative integer, got {}", v))?,
        None => 0,
    };
    let ids = match matches.value_of("calibration-ids") {
        Some(f) => Some(calibrate::read_ids(f)?),
        None => None,
    };
    let needed = NeededSites::from_annotations(regions, false);
    let allele_numbers = AlleleNumbers::read_from_file(an_file, &needed)?;
    eprintln!(
        "[INFO] allele numbers for {} of {} needed sites read from {}",
        allele_numbers.len(),
        needed.len(),
        an_file
    );
    let options = calibrate::Options {
        n_min,
        ids,
        filter_for_id: id,
    };
    let mut units = calibrate::collect_units(
        regions,
        genome,
        papa,
        &allele_numbers,
        observed_file,
        &needed,
        &options,
    )?;
    let (mut calibration, bins) = calibrate::fit(&mut units, n_bins, n_min)?;
    calibration.notes.push(format!(
        "source\tgenovo {} --action calibrate: {} possible synonymous point mutations, {} observed, {} bins of allele numbers",
        VERSION,
        units.len(),
        units.iter().filter(|u| u.y).count(),
        bins.len()
    ));
    calibration
        .notes
        .push("note\tS_indel is not calibrated; it is set to S_snv".to_string());
    calibration.write_to_file(out)?;
    let diagnostics = match matches.value_of("calibration-diagnostics") {
        Some(f) => Some(f.to_string()),
        None if out != "-" => Some(format!("{}.diagnostics.tsv", out)),
        None => None,
    };
    if let Some(d) = diagnostics {
        calibrate::write_diagnostics(&d, &units, &bins, &calibration)?;
    }
    eprintln!(
        "[INFO] S = {:e}, h at {} allele numbers from {} to {}",
        calibration.s_snv,
        calibration.knots().len(),
        calibration.knots()[0].0,
        calibration.knots()[calibration.knots().len() - 1].0
    );
    Ok(())
}

/// The scaling of mutation probabilities for `enumerate` and `classify`:
/// `--scaling-factor`, or `--calibration`, with `--allele-numbers` and
/// `--min-an`. Read once and kept for the rest of the run.
fn site_scaling(
    matches: &clap::ArgMatches,
    regions: &Option<Vec<mutexpect::SeqAnnotation>>,
    include_intronic: bool,
    scaling_factor: f32,
) -> Result<&'static SiteScaling> {
    static CELL: std::sync::OnceLock<SiteScaling> = std::sync::OnceLock::new();
    if let Some(s) = CELL.get() {
        return Ok(s);
    }
    let mut calibration = match matches.value_of("calibration") {
        Some(file) => Calibration::read_from_file(file)?,
        None => Calibration::uniform(scaling_factor as f64),
    };
    if let Some(min_an) = matches.value_of("min-an") {
        calibration.n_min = min_an
            .parse()
            .map_err(|_| anyhow::anyhow!("--min-an must be a non-negative integer, got {}", min_an))?;
    }
    let allele_numbers = match matches.value_of("allele-numbers") {
        Some(file) => {
            let regions = require_initialization(regions, "--genomic-regions or --gff3")?;
            let needed = NeededSites::from_annotations(regions, include_intronic);
            let an = AlleleNumbers::read_from_file(file, &needed)?;
            eprintln!(
                "[INFO] allele numbers for {} of {} needed sites read from {}",
                an.len(),
                needed.len(),
                file
            );
            Some(an)
        }
        None => None,
    };
    let scaling = SiteScaling::new(calibration, allele_numbers)?;
    Ok(CELL.get_or_init(|| scaling))
}

fn main() -> Result<()> {
    let app = App::new("genovo")
        .version(VERSION)
        .author("Jörn Bethune")
        .about("Determine genes enriched with de-novo mutations")
        .after_help("If no --action is given, all actions are executed.\n\
                     Possible actions are: transform, enumerate, expect, sample, classify, compare\n\
                     \n\
                     --action calibrate fits --calibration for population data with --allele-numbers on the\n\
                     possible synonymous mutations and the --observed-mutations (see the README).\n\
                     \n\
                     --action possible_mutations is not part of the pipeline: it writes every possible\n\
                     mutation (chr, position, ref, alt, transcript, mutation type, sequence context) to\n\
                     STDOUT. It needs --genome and --gff3 or --genomic-regions; the context has the\n\
                     k-mer size of --point-mutation-probabilities if given (otherwise 9).\n\
                     \n\
                     Progress messages are written to STDERR." )
        .arg(Arg::with_name("action")
             .long("action")
             .value_name("ACTION")
             .help("Only run a specific step in the pipeline")
             .takes_value(true))

        //raw data arguments
        .arg(Arg::with_name("gff3")
             .long("gff3")
             .value_name("FILE")
             .help("gff3 gene annotations file")
             .takes_value(true))
        .arg(Arg::with_name("genome")
             .long("genome")
             .value_name("FILE")
             .help("A 2bit reference genome sequence file")
             .takes_value(true))
        .arg(Arg::with_name("point-mutation-probabilities")
             .long("point-mutation-probabilities")
             .value_name("FILE")
             .help("A pattern partition prediction point mutation probability table")
             .takes_value(true))
        .arg(Arg::with_name("indel-mutation-probabilities")
             .long("indel-mutation-probabilities")
             .value_name("FILE")
             .help("A pattern partition prediction indel mutation probability table")
             .takes_value(true))
        .arg(Arg::with_name("observed-mutations")
             .long("observed-mutations")
             .value_name("FILE")
             .help("A vcf-like file containing observed point mutations")
             .takes_value(true))

        // input/output file arguments
        .arg(Arg::with_name("regions-bed")
             .long("regions-bed")
             .value_name("FILE")
             .help("With --action transform: also write the positions where allele numbers are needed (exons ± 2 bases, or whole regions with --include-intronic) as a BED file, e.g. to restrict a genome-wide allele number file")
             .takes_value(true))
        .arg(Arg::with_name("genomic-regions")
             .long("genomic-regions")
             .value_name("FILE")
             .help("Locations of exons, CDS and their phases for each gene")
             .takes_value(true))
        .arg(Arg::with_name("possible-mutations")
             .long("possible-mutations")
             .value_name("FILE")
             .help("A list of all possible point mutations for each gene")
             .takes_value(true))
        .arg(Arg::with_name("classified-mutations")
             .long("classified-mutations")
             .value_name("FILE")
             .help("Observed, classified point mutations")
             .takes_value(true))
        .arg(Arg::with_name("expected-mutations")
             .long("expected-mutations")
             .value_name("FILE")
             .help("Expected number of point mutations per gene")
             .takes_value(true))
        .arg(Arg::with_name("sampled-mutations")
             .long("sampled-mutations")
             .value_name("FILE")
             .help("Sampled number of point mutations per gene")
             .takes_value(true))
        .arg(Arg::with_name("significant-mutations")
             .long("significant-mutations")
             .value_name("FILE")
             .help("Statistical test results for every gene")
             .default_value("-")
             .takes_value(true))

        // non-file args
        .arg(Arg::with_name("id")
             .long("--id")
             .value_name("NAME")
             .help("Only process a gene/transcript with the given ID (useful for parallel processing)")
             .takes_value(true))
        .arg(Arg::with_name("scaling-factor")
             .long("--scaling-factor")
             .value_name("NAME")
             .help("Scaling factor for all mutation probabilities")
             .default_value("1.0")
             .takes_value(true))
        .arg(Arg::with_name("allele-numbers")
             .long("--allele-numbers")
             .value_name("FILE")
             .help("Number of called chromosomes at each site (tab-separated 'chrom pos an' or 'chrom:pos an', 1-based, sorted). Sites that are not in the file are not callable and are left out of expected and observed counts")
             .takes_value(true))
        .arg(Arg::with_name("calibration")
             .long("--calibration")
             .value_name("FILE")
             .help("Scaling factors and allele-number dependence of the mutation probabilities (instead of --scaling-factor); needs --allele-numbers if it depends on allele numbers. Output of --action calibrate")
             .conflicts_with("scaling-factor")
             .takes_value(true))
        .arg(Arg::with_name("calibration-diagnostics")
             .long("--calibration-diagnostics")
             .value_name("FILE")
             .help("With --action calibrate: observed and expected counts with the fitted calibration, by allele number and by rate (default: <calibration>.diagnostics.tsv)")
             .takes_value(true))
        .arg(Arg::with_name("calibration-bins")
             .long("--calibration-bins")
             .value_name("NUMBER")
             .help("With --action calibrate: number of allele-number bins for h")
             .default_value("25")
             .takes_value(true))
        .arg(Arg::with_name("calibration-ids")
             .long("--calibration-ids")
             .value_name("FILE")
             .help("With --action calibrate: only use the transcripts listed in this file (one ID per line), e.g. canonical transcripts")
             .takes_value(true))
        .arg(Arg::with_name("min-an")
             .long("--min-an")
             .value_name("NUMBER")
             .help("With --allele-numbers: leave out sites with a smaller allele number from expected and observed counts (overrides the minimum in --calibration)")
             .requires("allele-numbers")
             .takes_value(true))
        .arg(Arg::with_name("sum-up-observed-mutations-per-transcript")
             .long("--sum-up-observed-mutations-per-transcript")
             .help("Tally up the number of observed mutations instead of displaying them in detail (can be used in combination with --action classify and --classified-mutations <output-file>)"))
        .arg(Arg::with_name("include-intronic")
             .long("--include-intronic")
             .help("Also enumerate intronic variants. Default is to not enumerate intronic variants."))
        .arg(Arg::with_name("include-unknown")
             .long("--include-unknown")
             .help("Also enumerate variants with unknown type (primarily UTR). Default is to not enumerate variants with unknown type."))
        .arg(Arg::with_name("exclude-synonymous")
             .long("--exclude-synonymous")
             .help("Do not enumerate synonymous variants. Default is to enumerate synonymous variants."))
        .arg(Arg::with_name("number-of-random-samples")
             .long("--number-of-random-samples")
             .value_name("NUMBER")
             .default_value("1000")
             .help("The number of random samples that should be generated")
             .takes_value(true))
        .arg(Arg::with_name("seed")
             .long("--seed")
             .value_name("NUMBER")
             .help("Seed for the random numbers of the sample and compare steps. Runs with the same seed and inputs give identical results. Default: a random seed, which is reported on STDERR")
             .takes_value(true))
        .arg(Arg::with_name("required-tag")
             .long("--required-tag")
             .value_name("TAG")
             .help("Exclude entrences from gff file that to no include the tag TAG")
             .takes_value(true)
             .multiple(true))
        .arg(Arg::with_name("filter-plof")
             .long("--filter-plof")
             .help("filter putative LoF variants. (50 bp rule)"))
        .arg(Arg::with_name("alpha")
             .long("--alpha")
             .value_name("FLOAT")
             .default_value("0.1")
             .help("Alpha value used for confidence interval. 0.1 -> 90% c.i.")
             .takes_value(true))
    ;


    let matches = app.get_matches();
    let run_all = matches.value_of("action").is_none();
    let id = matches.value_of("id");
    let scaling_factor: f32 = matches
        .value_of("scaling-factor")
        .expect("default value")
        .parse()?;
    let include_intronic = matches.occurrences_of("include-intronic") > 0;
    let include_unknown = matches.occurrences_of("include-unknown") > 0;
    let filter_plof = matches.occurrences_of("filter-plof") > 0;
  
    let seed: u64 = match matches.value_of("seed") {
        Some(s) => s
            .parse()
            .map_err(|_| anyhow::anyhow!("--seed must be a non-negative integer, got {}", s))?,
        None => {
            let seed = rand::random();
            if run_all || matches!(matches.value_of("action"), Some("sample") | Some("compare")) {
                eprintln!("Random seed: {} (use --seed {} to reproduce this run)", seed, seed);
            }
            seed
        }
    };

    let required_tags: Option<Vec<&str>> = {
       if let Some(tags) = matches.values_of("required-tag"){
           Some(tags.collect())
       } else {
           None
       }
    };
   
    //println!("{:#?}", required_tags);
    
    
    /*
     * Depending on what --action is given, each step may or may not produce results.
     * Therefore the variables are all Option's.
     */

    let ref_genome = {
        if let Some(ref_genome_file) = matches.value_of("genome") {
            Some(TwoBitFile::open(ref_genome_file, false)?)
        } else {
            None
        }
    };

    let papa = {
        if let Some(papa_file) = matches.value_of("point-mutation-probabilities") {
            Some(PaPaPred::new(papa_file, Some(5))?) // 5 to have at least 2 flanking bases around a point mutation to always have a full codon available for every coding site
        } else {
            None
        }
    };

    let papa_indel = {
        if let Some(papa_file) = matches.value_of("indel-mutation-probabilities") {
            let min_width = {
                if let Some(papa_point) = &papa {
                    Some(papa_point.kmer_size().saturating_sub(1))
                } else {
                    None
                }
            };
            Some(PaPaPredIndel::new(papa_file, min_width)?)
        } else {
            None
        }
    };

    // (only classify reads the observed mutations into memory; calibrate
    // streams them)
    let observed_mutations = {
        if !(run_all || matches.value_of("action") == Some("classify")) {
            None
        } else if let Some(observed_mutations_file) = matches.value_of("observed-mutations") {
            Some(read_observed_mutations_from_file(
                observed_mutations_file,
                -1,
            )?) //TODO expose adjustment parameter to CLI
        } else {
            None
        }
    };

    eprintln!("transform");
    // action=transform
    let regions = {
        if let Some(gff3) = matches.value_of("gff3") {
            let regions = transform::transform_gff3_annotations(gff3, id, required_tags)?;
            if let Some(regions_file) = matches.value_of("genomic-regions") {
                transform::write_sequence_annotations_to_file(regions_file, &regions)?;
            }
            if let Some(bed) = matches.value_of("regions-bed") {
                NeededSites::from_annotations(&regions, include_intronic).write_bed(bed)?;
            }
            if matches.value_of("action") == Some("transform") {
                // we are done here
                return Ok(());
            }
            Some(regions)
        } else if let Some(regions_file) = matches.value_of("genomic-regions") {
            let regions = mutexpect::read_sequence_annotations_from_file(regions_file, id)?;
            if matches.value_of("action") == Some("transform") {
                // only the BED file is left to do
                let bed = matches.value_of("regions-bed").ok_or_else(|| {
                    anyhow::anyhow!("--action transform needs --gff3 (or --genomic-regions with --regions-bed)")
                })?;
                NeededSites::from_annotations(&regions, include_intronic).write_bed(bed)?;
                return Ok(());
            }
            Some(regions)
        } else if run_all || matches.value_of("action") == Some("transform") {
            {
                return Err(anyhow::anyhow!("Please provide the --gff3 parameter"));
            }
        } else {
            None
        }
    };

    
    //action=calibrate
    if matches.value_of("action") == Some("calibrate") {
        eprintln!("calibrate");
        run_calibrate(&matches, &regions, &ref_genome, &papa, id)?;
        return Ok(());
    }

    //action=enumerate
    let possible_mutations = {
        if run_all || matches.value_of("action") == Some("enumerate") {
            eprintln!("enumerate");
            let possible_mutations = enumerate_possible_mutations(
                require_initialization(&regions, "--genomic-regions")?,
                require_initialization(&ref_genome, "--genome")?,
                require_initialization(&papa, "--point-mutation-probabilities")?,
                &papa_indel,
                site_scaling(&matches, &regions, include_intronic, scaling_factor)?,
                true,
                id,
                include_intronic,
                include_unknown,
                filter_plof,
            )?;

            if let Some(possible_mutations_file) = matches.value_of("possible-mutations") {
                enumerate::write_to_file(possible_mutations_file, &possible_mutations)?
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(possible_mutations)
        } else if let Some(possible_mutations_file) = matches.value_of("possible-mutations") {
            Some(enumerate::read_from_file(possible_mutations_file)?)
        } else {
            None
        }
    };


    //action=possible_mutations
    if matches.value_of("action") == Some("possible_mutations") {
        eprintln!("possible_mutations");
        let radius = papa.as_ref().map_or(4, |p| p.radius());
        let stdout = std::io::stdout();
        write_possible_mutations(
            std::io::BufWriter::new(stdout.lock()),
            require_initialization(&regions, "--genomic-regions or --gff3")?,
            require_initialization(&ref_genome, "--genome")?,
            radius,
            id,
            include_intronic,
            include_unknown,
            filter_plof,
        )?;
        return Ok(());
    };

    eprintln!("expect");
    //action=expect
    let expected_mutations = {
        if run_all || matches.value_of("action") == Some("expect") {
            let expected_mutations = expected_number_of_mutations(
                require_initialization(&possible_mutations, "--possible-mutations")?,
                id,
            )?;
            if let Some(expected_mutations_file) = matches.value_of("expected-mutations") {
                expect::write_to_file(expected_mutations_file, &expected_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(expected_mutations)
        } else if let Some(expected_mutations_file) = matches.value_of("expected-mutations") {
            Some(expect::read_from_file(expected_mutations_file)?)
        } else {
            None
        }
    };

    eprintln!("sample");
    //action=sample
    let sampled_mutations = {
        if run_all || matches.value_of("action") == Some("sample") {
            let iterations: usize = matches
                .value_of("number-of-random-samples")
                .expect("clap default value")
                .parse()
                .unwrap(); //TODO proper error handling
            let sampled_mutations = sample_mutations(
                require_initialization(&possible_mutations, "--possible-mutations")?,
                iterations,
                id,
                seed,
            )?;

            if let Some(sampled_mutations_file) = matches.value_of("sampled-mutations") {
                sample::write_to_file(sampled_mutations_file, &sampled_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(sampled_mutations)
        } else if let Some(sampled_mutations_file) = matches.value_of("sampled-mutations") {
            Some(sample::read_from_file(sampled_mutations_file)?)
        } else {
            None
        }
    };

    std::mem::drop(possible_mutations); // let's free up some memory

    eprintln!("classify");
    //action=classify
    let classified_mutations = {
        if run_all || matches.value_of("action") == Some("classify") {
            let classified_mutations = classify_mutations(
                require_initialization(&observed_mutations, "--observed-mutations")?,
                require_initialization(&regions, "--genomic-regions")?,
                require_initialization(&ref_genome, "--genome")?,
                id,
                filter_plof,
                site_scaling(&matches, &regions, include_intronic, scaling_factor)?,
            )?;

            if let Some(classified_mutations_file) = matches.value_of("classified-mutations") {
                if matches.occurrences_of("sum-up-observed-mutations-per-transcript") > 0 {
                    // sum up the number of mutations of each type for each transcript
                    observed::sum_up_and_write_to_file(
                        classified_mutations_file,
                        &classified_mutations,
                    )?;
                } else {
                    // Write the classification for each individual mutation on a separate line
                    // This is what you want if you want to run the full pipeline
                    observed::write_to_file(classified_mutations_file, &classified_mutations)?;
                }
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(classified_mutations)
        } else if let Some(classified_mutations_file) = matches.value_of("classified-mutations") {
            Some(observed::read_from_file(classified_mutations_file)?)
        } else {
            None
        }
    };

    eprintln!("compare");
    //action=compare
    let _significant_mutations = {
        if run_all || matches.value_of("action") == Some("compare") {
            let alpha: Float = matches
                .value_of("alpha")
                .expect("clap default value")
                .parse()
                .unwrap(); //TODO proper error handling
            let significant_mutations = compare_mutations(
                require_initialization(&classified_mutations, "--classified-mutations")?,
                require_initialization(&expected_mutations, "--expected-mutations")?,
                require_initialization(&sampled_mutations, "--sampled-mutations")?,
                id,
                alpha,
                seed,
            )?;

            if let Some(significant_mutations_file) = matches.value_of("significant-mutations") {
                compare::write_to_file(significant_mutations_file, &significant_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(significant_mutations)
        } else {
            // last step in the pipeline. No work needs to be done
            None
        }
    };

    if !run_all {
        return Err(anyhow::anyhow!(
            "Invalid --action parameter: {}",
            matches.value_of("action").expect("clap")
        ));
    }
    Ok(())
}
