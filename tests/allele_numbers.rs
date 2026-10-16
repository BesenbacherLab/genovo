//! `--allele-numbers`, `--calibration` and `--min-an` on the chr22 test data
//! (see `docs/design/allele-number-scaling.md`), for one transcript with
//! observed LoF mutations, using synthetic allele numbers.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::Once;

const ID: &str = "ENST00000263253.7";
const GFF3: &str = "gencode.v19.CCDS.annotation.protein_coding_chr22.gff3.gz";
const GENOME: &str = "hg19_chr22.2bit";
const OBSERVED: &str = "observed_mutations_DDD_2017_chr22.txt";
const SCALING: &str = "8586";

fn test_data() -> PathBuf {
    static UNPACK: Once = Once::new();
    // not shared with tests/pipeline.rs, which may unpack at the same time
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genovo_test_data_an");
    UNPACK.call_once(|| {
        if !dir.join(GENOME).exists() {
            fs::create_dir_all(&dir).unwrap();
            let tarball = Path::new(env!("CARGO_MANIFEST_DIR")).join("test_data/test_data.tar.gz");
            let ok = Command::new("tar")
                .arg("xzf")
                .arg(&tarball)
                .arg("-C")
                .arg(&dir)
                .status()
                .expect("run tar")
                .success();
            assert!(ok, "unpacking {}", tarball.display());
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

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_genovo"))
        .current_dir(test_data())
        .args(args)
        .output()
        .expect("run genovo")
}

fn genovo(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "genovo {:?} failed:\n{}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stderr).unwrap()
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

fn model() -> Vec<&'static str> {
    vec![
        "--point-mutation-probabilities",
        "SNV_mutation_rate_model.txt",
        "--indel-mutation-probabilities",
        "indel_mutation_rate_model.txt",
    ]
}

struct Outputs {
    results: String,
    expected: String,
    sampled: String,
    classified: String,
}

/// The whole pipeline for `ID` in one go.
fn run_all(dir: &Path, extra: &[&str]) -> Outputs {
    let f = |n: &str| dir.join(n);
    let mut args = vec![
        "--gff3", GFF3, "--genome", GENOME, "--observed-mutations", OBSERVED, "--id", ID,
        "--seed", "1",
    ];
    let (res, exp, sam, pos, cls) = (
        f("results.tsv"),
        f("expected.tsv"),
        f("sampled.tsv"),
        f("possible.txt"),
        f("classified.tsv"),
    );
    args.extend(&[
        "--significant-mutations", p(&res), "--expected-mutations", p(&exp),
        "--sampled-mutations", p(&sam), "--possible-mutations", p(&pos),
        "--classified-mutations", p(&cls),
    ]);
    args.extend(model());
    args.extend(extra);
    genovo(&args);
    Outputs {
        results: fs::read_to_string(res).unwrap(),
        expected: fs::read_to_string(exp).unwrap(),
        sampled: fs::read_to_string(sam).unwrap(),
        classified: fs::read_to_string(cls).unwrap(),
    }
}

/// enumerate, expect, sample, classify and compare as separate runs.
fn run_steps(dir: &Path, extra: &[&str]) -> Outputs {
    let f = |n: &str| dir.join(n).to_str().unwrap().to_string();
    let (pos, exp, sam, cls, res) = (
        f("possible.txt"),
        f("expected.tsv"),
        f("sampled.tsv"),
        f("classified.tsv"),
        f("results.tsv"),
    );
    let step = |args: &[&str]| {
        let mut a = vec!["--id", ID];
        a.extend(args);
        a.extend(extra);
        genovo(&a);
    };
    let mut enumerate = vec![
        "--action", "enumerate", "--gff3", GFF3, "--genome", GENOME, "--possible-mutations", &pos,
    ];
    enumerate.extend(model());
    step(&enumerate);
    step(&["--action", "expect", "--possible-mutations", &pos, "--expected-mutations", &exp]);
    step(&[
        "--action", "sample", "--possible-mutations", &pos, "--sampled-mutations", &sam, "--seed",
        "1",
    ]);
    step(&[
        "--action", "classify", "--gff3", GFF3, "--genome", GENOME, "--observed-mutations",
        OBSERVED, "--classified-mutations", &cls,
    ]);
    step(&[
        "--action", "compare", "--classified-mutations", &cls, "--expected-mutations", &exp,
        "--sampled-mutations", &sam, "--significant-mutations", &res, "--seed", "1",
    ]);
    Outputs {
        results: fs::read_to_string(res).unwrap(),
        expected: fs::read_to_string(exp).unwrap(),
        sampled: fs::read_to_string(sam).unwrap(),
        classified: fs::read_to_string(cls).unwrap(),
    }
}

fn stored(kind: &str) -> String {
    let f = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(format!("{}.{}.tsv", ID, kind));
    fs::read_to_string(f).unwrap()
}

/// 0-based positions where allele numbers are needed for `ID`.
fn needed_positions(dir: &Path) -> Vec<(String, usize)> {
    let bed = dir.join("needed.bed");
    genovo(&["--action", "transform", "--gff3", GFF3, "--id", ID, "--regions-bed", p(&bed)]);
    let mut out = Vec::new();
    for line in fs::read_to_string(bed).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (start, stop): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        out.extend((start..stop).map(|pos| (f[0].to_string(), pos)));
    }
    out
}

/// Write an allele number file (gnomAD `chrom:pos<TAB>AN` layout, 1-based)
/// with `an(pos)` at the needed 0-based positions; `None` leaves a site out.
fn write_allele_numbers(dir: &Path, name: &str, an: impl Fn(usize) -> Option<u32>) -> PathBuf {
    let mut text = String::from("locus\tAN\n");
    for (chr, pos) in needed_positions(dir) {
        if let Some(n) = an(pos) {
            text.push_str(&format!("{}:{}\t{}\n", chr, pos + 1, n));
        }
    }
    let path = dir.join(name);
    fs::write(&path, text).unwrap();
    path
}

/// 0-based allele-number sites of the observed mutations of `ID`: the
/// mutated base, or the base after the anchor for indels.
fn observed_sites(classified: &str) -> Vec<usize> {
    classified
        .lines()
        .skip(1)
        .filter(|l| l.starts_with(ID))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let pos: usize = f[2].parse().unwrap();
            let (from, to) = f[4].split_once("->").unwrap();
            if from.len() == 1 && to.len() == 1 {
                pos
            } else {
                pos + 1
            }
        })
        .collect()
}

/// Expected count per (region, mutation type).
fn expected_counts(results: &str) -> Vec<(String, f64, usize)> {
    results
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[1].to_string(), f[3].parse().unwrap(), f[2].parse().unwrap())
        })
        .collect()
}

#[test]
fn allele_numbers_without_calibration_change_nothing_if_all_sites_are_callable() {
    let dir = out_dir("genovo_an_uniform");
    let an = write_allele_numbers(&dir, "an.tsv", |_| Some(100_000));
    let o = run_all(&dir, &["--scaling-factor", SCALING, "--allele-numbers", p(&an)]);
    assert_eq!(o.results, stored("results"));
    assert_eq!(o.expected, stored("expected"));
    assert_eq!(o.sampled, stored("sampled"));
}

#[test]
fn calibration_h_scales_the_exponent() {
    // S = 2 × 8586 and h = 0.5 at every site: the same exponent as
    // --scaling-factor 8586, so identical results.
    let dir = out_dir("genovo_an_calibration");
    let an = write_allele_numbers(&dir, "an.tsv", |_| Some(1000));
    let calibration = dir.join("calibration.tsv");
    fs::write(
        &calibration,
        "# genovo calibration v1\n# S_snv\t17172\n# S_indel\t17172\n# n_min\t0\nan\th\n1000\t0.5\n",
    )
    .unwrap();
    let o = run_all(&dir, &["--calibration", p(&calibration), "--allele-numbers", p(&an)]);
    assert_eq!(o.results, stored("results"));
    assert_eq!(o.expected, stored("expected"));
    assert_eq!(o.sampled, stored("sampled"));
    // Without allele numbers, this calibration is refused.
    let mut args = vec!["--action", "enumerate", "--gff3", GFF3, "--genome", GENOME, "--id", ID];
    args.extend(model());
    args.extend(&["--calibration", p(&calibration)]);
    let out = run(&args);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--allele-numbers"));
}

#[test]
fn masked_sites_leave_expected_and_observed_counts() {
    let dir = out_dir("genovo_an_masked");
    let unmasked = run_all(&out_dir("genovo_an_unmasked"), &["--scaling-factor", SCALING]);
    let sites = observed_sites(&unmasked.classified);
    assert!(sites.len() >= 12, "{:?}", sites);

    // The observed sites get a low allele number and are masked with --min-an.
    let low = write_allele_numbers(&dir, "low.tsv", |pos| {
        Some(if sites.contains(&pos) { 10 } else { 100_000 })
    });
    let args = ["--scaling-factor", SCALING, "--allele-numbers", p(&low), "--min-an", "100"];
    let masked = run_all(&dir, &args);
    assert!(observed_sites(&masked.classified).is_empty(), "{}", masked.classified);
    let before = expected_counts(&unmasked.results);
    let after = expected_counts(&masked.results);
    assert_eq!(before.len(), after.len());
    for (t, e1, o1) in &after {
        assert_eq!(*o1, 0, "{}: observed mutations at masked sites", t);
        let e0 = before.iter().find(|(t0, _, _)| t0 == t).unwrap().1;
        assert!(*e1 <= e0 * (1.0 + 1e-6), "{}: {} > {}", t, e1, e0);
    }
    let total = |c: &[(String, f64, usize)]| c.iter().map(|(_, e, _)| e).sum::<f64>();
    assert!(total(&after) < total(&before));

    // Leaving the sites out of the file (not callable) masks them too.
    let missing = write_allele_numbers(&dir, "missing.tsv", |pos| {
        if sites.contains(&pos) {
            None
        } else {
            Some(100_000)
        }
    });
    let dir2 = out_dir("genovo_an_missing");
    let o = run_all(&dir2, &["--scaling-factor", SCALING, "--allele-numbers", p(&missing)]);
    assert_eq!(o.expected, masked.expected);
    assert_eq!(o.sampled, masked.sampled);
    assert_eq!(o.results, masked.results);

    // Step by step gives the same result.
    let steps = run_steps(&out_dir("genovo_an_masked_steps"), &args);
    assert_eq!(steps.results, masked.results);
    assert_eq!(steps.expected, masked.expected);
    assert_eq!(steps.sampled, masked.sampled);
    assert_eq!(steps.classified, masked.classified);
}

#[test]
fn option_errors() {
    let dir = out_dir("genovo_an_errors");
    let calibration = dir.join("c.tsv");
    fs::write(&calibration, "# genovo calibration v1\n# S_snv\t1\nan\th\n").unwrap();
    let base = ["--action", "enumerate", "--gff3", GFF3, "--genome", GENOME, "--id", ID];
    for extra in [
        vec!["--min-an", "10"],                                            // needs --allele-numbers
        vec!["--calibration", p(&calibration), "--scaling-factor", "2"],   // conflict
        vec!["--allele-numbers", "does-not-exist.tsv"],
    ] {
        let mut args = base.to_vec();
        args.extend(model());
        args.extend(&extra);
        assert!(!run(&args).status.success(), "{:?} should fail", extra);
    }
}
