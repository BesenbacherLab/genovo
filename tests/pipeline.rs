//! End-to-end tests on the chr22 test data set (`test_data/test_data.tar.gz`).
//!
//! * The classification of observed mutations must agree with the
//!   classification of the same mutations when they are enumerated as possible
//!   mutations (otherwise observed and expected counts measure different things).
//! * With `--seed`, the pipeline output is compared with stored results
//!   (`tests/data/`), run in one go and step by step.
//!
//! To update the stored results after an intended change, run the tests with
//! `GENOVO_BLESS=1` and review the diff of `tests/data/`. The full chr22 test
//! is ignored by default (slow in debug builds):
//! `cargo test --release --test pipeline -- --ignored`.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Once;

/// Transcripts used by the fast tests: plus strand with observed LoF
/// mutations, and two minus-strand transcripts (one with 11 exons).
const IDS: [&str; 3] = ["ENST00000263253.7", "ENST00000215862.4", "ENST00000343518.6"];

const SEED: &str = "1";

fn test_data() -> PathBuf {
    static UNPACK: Once = Once::new();
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genovo_test_data");
    UNPACK.call_once(|| {
        if !dir.join("hg19_chr22.2bit").exists() {
            fs::create_dir_all(&dir).unwrap();
            let tarball = Path::new(env!("CARGO_MANIFEST_DIR")).join("test_data/test_data.tar.gz");
            let status = Command::new("tar")
                .arg("xzf")
                .arg(&tarball)
                .arg("-C")
                .arg(&dir)
                .status()
                .expect("run tar");
            assert!(status.success(), "unpacking {}", tarball.display());
        }
    });
    dir
}

fn out_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

/// Run genovo in the test data directory; returns stdout.
fn genovo(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_genovo"))
        .current_dir(test_data())
        .args(args)
        .output()
        .expect("run genovo");
    assert!(
        out.status.success(),
        "genovo {:?} failed:\n{}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn path(p: &Path) -> &str {
    p.to_str().unwrap()
}

const GFF3: &str = "gencode.v19.CCDS.annotation.protein_coding_chr22.gff3.gz";
const GENOME: &str = "hg19_chr22.2bit";
const OBSERVED: &str = "observed_mutations_DDD_2017_chr22.txt";

fn model_args() -> Vec<&'static str> {
    vec![
        "--point-mutation-probabilities",
        "SNV_mutation_rate_model.txt",
        "--indel-mutation-probabilities",
        "indel_mutation_rate_model.txt",
        "--scaling-factor",
        "8586",
    ]
}

/// Mutation type names as written by `classify` (numeric codes).
fn type_name(code: &str) -> &'static str {
    match code {
        "0" => "Unknown",
        "1" => "Synonymous",
        "2" => "Missense",
        "3" => "Nonsense",
        "4" => "StopLoss",
        "5" => "StartCodon",
        "6" => "SpliceSite",
        "7" => "Intronic",
        "8" => "InFrameIndel",
        "9" => "FrameshiftIndel",
        _ => panic!("unknown mutation type code {}", code),
    }
}

/// For each transcript in `IDS`: (expected type, observed type) for every
/// possible mutation whose two classifications differ, and the number compared.
fn classification_differences(scratch: &str) -> (BTreeMap<(String, String), usize>, usize) {
    let dir = out_dir(scratch);
    let mut differences = BTreeMap::new();
    let mut compared = 0;
    for id in IDS {
        // All possible mutations, then the same mutations as "observed" ones.
        let listing = genovo(&[
            "--action",
            "possible_mutations",
            "--gff3",
            GFF3,
            "--genome",
            GENOME,
            "--point-mutation-probabilities",
            "SNV_mutation_rate_model.txt",
            "--id",
            id,
        ]);
        let mut expected = BTreeMap::new();
        let mut observed_input = String::new();
        for line in listing.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            let (chr, pos, reference, alt, region, mtype) = (f[0], f[1], f[2], f[3], f[4], f[5]);
            assert_eq!(region, id);
            // a concrete indel of the right kind: in-frame (3) or frameshift (1)
            let alt = alt.replace('N', "A");
            observed_input.push_str(&format!("{} {} {} {}\n", chr, pos, reference, alt));
            expected.insert((pos.to_string(), alt), mtype.to_string());
        }
        let observed_file = dir.join(format!("{}.observed.txt", id));
        fs::write(&observed_file, observed_input).unwrap();
        let classified_file = dir.join(format!("{}.classified.tsv", id));
        genovo(&[
            "--action",
            "classify",
            "--gff3",
            GFF3,
            "--genome",
            GENOME,
            "--observed-mutations",
            path(&observed_file),
            "--classified-mutations",
            path(&classified_file),
            "--id",
            id,
        ]);
        let classified = fs::read_to_string(&classified_file).unwrap();
        let mut observed = BTreeMap::new();
        for line in classified.lines().skip(1) {
            // region, chromosome, 0-based position, type code, ref->alt
            let f: Vec<&str> = line.split('\t').collect();
            let pos = (f[2].parse::<usize>().unwrap() + 1).to_string();
            let alt = f[4].split("->").nth(1).unwrap().to_string();
            observed.insert((pos, alt), type_name(f[3]).to_string());
        }
        assert_eq!(observed.len(), expected.len(), "{}", id);
        for (key, e) in &expected {
            let o = observed.get(key).map(String::as_str).unwrap_or("missing");
            compared += 1;
            if o != e {
                *differences.entry((e.clone(), o.to_string())).or_insert(0) += 1;
            }
        }
    }
    (differences, compared)
}

/// Every possible mutation, given to `classify` as an observed mutation,
/// gets the type `enumerate` gave it.
#[test]
fn observed_classification_matches_possible_mutations() {
    let (differences, compared) = classification_differences("genovo_consistency");
    assert!(compared > 10_000, "only {} mutations compared", compared);
    assert!(
        differences.is_empty(),
        "(enumerate, classify) -> number of possible mutations: {:?}",
        differences
    );
}

/// Compare `actual` with the stored file `tests/data/<name>`, or store it
/// with `GENOVO_BLESS=1`.
fn check_golden(name: &str, actual: &str) {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data").join(name);
    if std::env::var_os("GENOVO_BLESS").is_some() {
        fs::create_dir_all(golden.parent().unwrap()).unwrap();
        fs::write(&golden, actual).unwrap();
        return;
    }
    let want = fs::read_to_string(&golden).unwrap_or_else(|_| {
        panic!(
            "{} is missing; run with GENOVO_BLESS=1 to create it",
            golden.display()
        )
    });
    if want != actual {
        let first = want
            .lines()
            .zip(actual.lines())
            .position(|(a, b)| a != b)
            .unwrap_or_else(|| want.lines().count().min(actual.lines().count()));
        panic!(
            "{} differs from the stored result (first difference at line {}):\n  stored: {:?}\n  now:    {:?}\n\
             If the change is intended, rerun with GENOVO_BLESS=1 and review the diff.",
            name,
            first + 1,
            want.lines().nth(first),
            actual.lines().nth(first)
        );
    }
}

/// Run the whole pipeline in one go; returns (results, expected, sampled).
fn run_all(dir: &Path, extra: &[&str]) -> (String, String, String) {
    let (res, exp, sam, pos) = (
        dir.join("results.tsv"),
        dir.join("expected.tsv"),
        dir.join("sampled.tsv"),
        dir.join("possible.txt"),
    );
    let mut args = vec![
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--observed-mutations",
        OBSERVED,
        "--seed",
        SEED,
        "--significant-mutations",
        path(&res),
        "--expected-mutations",
        path(&exp),
        "--sampled-mutations",
        path(&sam),
        "--possible-mutations",
        path(&pos),
    ];
    args.extend(model_args());
    args.extend(extra);
    genovo(&args);
    (
        fs::read_to_string(res).unwrap(),
        fs::read_to_string(exp).unwrap(),
        fs::read_to_string(sam).unwrap(),
    )
}

/// Run the pipeline one step at a time; returns (results, expected, sampled).
fn run_steps(dir: &Path, extra: &[&str]) -> (String, String, String) {
    let p = |name: &str| dir.join(name).to_str().unwrap().to_string();
    let (regions, possible, expected, sampled, classified, results) = (
        p("regions.txt"),
        p("possible.txt"),
        p("expected.tsv"),
        p("sampled.tsv"),
        p("classified.tsv"),
        p("results.tsv"),
    );
    let step = |args: &[&str]| {
        let mut a = args.to_vec();
        a.extend(extra);
        genovo(&a);
    };
    step(&["--action", "transform", "--gff3", GFF3, "--genomic-regions", &regions]);
    let mut enumerate = vec![
        "--action",
        "enumerate",
        "--genomic-regions",
        &regions,
        "--genome",
        GENOME,
        "--possible-mutations",
        &possible,
    ];
    enumerate.extend(model_args());
    step(&enumerate);
    step(&["--action", "expect", "--possible-mutations", &possible, "--expected-mutations", &expected]);
    step(&[
        "--action", "sample", "--possible-mutations", &possible, "--sampled-mutations", &sampled,
        "--seed", SEED,
    ]);
    step(&[
        "--action", "classify", "--genomic-regions", &regions, "--genome", GENOME,
        "--observed-mutations", OBSERVED, "--classified-mutations", &classified,
    ]);
    step(&[
        "--action", "compare", "--classified-mutations", &classified, "--expected-mutations",
        &expected, "--sampled-mutations", &sampled, "--significant-mutations", &results,
        "--seed", SEED,
    ]);
    (
        fs::read_to_string(results).unwrap(),
        fs::read_to_string(expected).unwrap(),
        fs::read_to_string(sampled).unwrap(),
    )
}

#[test]
fn pipeline_selected_transcripts() {
    for id in IDS {
        let one_go = run_all(&out_dir(&format!("genovo_all_{}", id)), &["--id", id]);
        let steps = run_steps(&out_dir(&format!("genovo_steps_{}", id)), &["--id", id]);
        assert!(one_go == steps, "{}: step-by-step run differs from the run in one go", id);
        check_golden(&format!("{}.results.tsv", id), &one_go.0);
        check_golden(&format!("{}.expected.tsv", id), &one_go.1);
        check_golden(&format!("{}.sampled.tsv", id), &one_go.2);
    }
}

#[test]
#[ignore = "slow in debug builds; run with cargo test --release --test pipeline -- --ignored"]
fn pipeline_chr22() {
    let one_go = run_all(&out_dir("genovo_all_chr22"), &[]);
    check_golden("chr22.results.tsv", &one_go.0);
    check_golden("chr22.expected.tsv", &one_go.1);
    check_golden("chr22.sampled.tsv", &one_go.2);
    let steps = run_steps(&out_dir("genovo_steps_chr22"), &[]);
    assert!(one_go == steps, "step-by-step run differs from the run in one go");
    // A single transcript gives the same rows as the whole chromosome.
    for id in IDS {
        let single = run_all(&out_dir(&format!("genovo_single_{}", id)), &["--id", id]);
        let rows = |t: &str| -> Vec<String> {
            let mut v: Vec<String> = t.lines().filter(|l| l.starts_with(id)).map(String::from).collect();
            v.sort();
            v
        };
        assert_eq!(rows(&single.0), rows(&one_go.0), "{} results", id);
        assert_eq!(rows(&single.1), rows(&one_go.1), "{} expected", id);
        assert_eq!(rows(&single.2), rows(&one_go.2), "{} sampled", id);
    }
}

/// A transcript whose sequence contexts reach into an N gap (chr22 starts
/// with Ns up to 0-based position 16050000 in the test genome) is skipped
/// with a warning by enumerate and by possible_mutations with probabilities,
/// and the other transcripts are unaffected.
#[test]
fn transcripts_next_to_n_are_skipped() {
    let dir = out_dir("genovo_n_gap");
    let regions = dir.join("regions.txt");
    let real = dir.join("real.txt");
    genovo(&["--action", "transform", "--gff3", GFF3, "--id", IDS[0], "--genomic-regions", path(&real)]);
    let real_line = fs::read_to_string(&real).unwrap();
    fs::write(
        &regions,
        format!("NEAR_GAP\tchr22\t+\t16050000-16050200\t16050000-16050200\t16050003-16050198\t0\n{}", real_line),
    )
    .unwrap();
    let expected = dir.join("expected.tsv");
    let possible = dir.join("possible.txt");
    let mut args = vec!["--action", "enumerate", "--genomic-regions", path(&regions), "--genome", GENOME,
                        "--possible-mutations", path(&possible)];
    args.extend(model_args());
    genovo(&args);
    genovo(&["--action", "expect", "--possible-mutations", path(&possible), "--expected-mutations", path(&expected)]);
    let expected = fs::read_to_string(&expected).unwrap();
    assert!(!expected.contains("NEAR_GAP"));
    assert!(expected.contains(IDS[0]));

    let mut args = vec!["--action", "possible_mutations", "--genomic-regions", path(&regions), "--genome", GENOME];
    args.extend(model_args());
    let sites = genovo(&args);
    assert!(!sites.contains("NEAR_GAP"));
    let total: f64 = sites.lines().map(|l| l.split('\t').nth(7).unwrap().parse::<f64>().unwrap()).sum();
    let mut lines = expected.lines();
    lines.next();
    let expected_total: f64 = lines
        .flat_map(|l| l.split('\t').skip(1).map(|v| v.parse::<f64>().unwrap()).collect::<Vec<_>>())
        .sum();
    assert!((total - expected_total).abs() < 1e-5 * expected_total, "{} vs {}", total, expected_total);
}
