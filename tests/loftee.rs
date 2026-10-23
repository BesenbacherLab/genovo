//! `--loftee`, `--loftee-filter` and `--gerp` on the chr22 test data, with
//! synthetic GERP++ scores.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::Once;

/// Plus strand with observed LoF mutations, and two minus-strand transcripts.
const IDS: [&str; 3] = ["ENST00000263253.7", "ENST00000215862.4", "ENST00000343518.6"];
const GFF3: &str = "gencode.v19.CCDS.annotation.protein_coding_chr22.gff3.gz";
const GENOME: &str = "hg19_chr22.2bit";
const OBSERVED: &str = "observed_mutations_DDD_2017_chr22.txt";

fn test_data() -> PathBuf {
    static UNPACK: Once = Once::new();
    // not shared with the other test files, which may unpack at the same time
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("genovo_test_data_loftee");
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

/// Run genovo, which must succeed; returns STDOUT.
fn stdout(args: &[&str]) -> String {
    let out = run(args);
    assert!(
        out.status.success(),
        "genovo {:?} failed:\n{}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

/// A GERP++ score file (`chrom pos score`, 1-based) for all transcripts'
/// exons: blocks of 7 bases scoring 5, then 14 scoring -2.
fn write_gerp(dir: &Path) -> PathBuf {
    let bed = dir.join("needed.bed");
    stdout(&["--action", "transform", "--gff3", GFF3, "--regions-bed", p(&bed)]);
    let mut text = String::from("chrom\tpos\tscore\n");
    for line in fs::read_to_string(bed).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let (start, stop): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
        for pos in start + 1..=stop {
            let score = if (pos / 7) % 3 == 0 { 5.0 } else { -2.0 };
            text.push_str(&format!("{}\t{}\t{}\n", f[0], pos, score));
        }
    }
    let path = dir.join("gerp.tsv");
    fs::write(&path, text).unwrap();
    path
}

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

/// Every possible mutation, given to `classify` as an observed mutation, gets
/// the LOFTEE label and numbers that `possible_mutations` gave it.
#[test]
fn labels_agree_between_possible_mutations_and_classify() {
    let dir = out_dir("genovo_loftee_consistency");
    let gerp = write_gerp(&dir);
    let mut labels: BTreeMap<String, usize> = BTreeMap::new();
    let mut fifty_bp_rule_fails_but_conserved = 0;
    for id in IDS {
        let loftee = ["--loftee", "--gerp", p(&gerp)];
        let mut args = vec![
            "--action",
            "possible_mutations",
            "--gff3",
            GFF3,
            "--genome",
            GENOME,
            "--id",
            id,
        ];
        args.extend(loftee);
        let listing = stdout(&args);
        let mut expected = BTreeMap::new();
        let mut observed_input = String::new();
        for line in listing.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 9, "{}", line);
            let (chr, pos, reference, alt, mtype, label, info) = (f[0], f[1], f[2], f[3], f[5], f[7], f[8]);
            let lof = matches!(mtype, "Nonsense" | "FrameshiftIndel" | "SpliceSite");
            assert_eq!(label != ".", lof, "{}", line);
            if lof {
                *labels.entry(format!("{} {}", mtype, label)).or_insert(0) += 1;
                if label == "HC" && info.contains("50_BP_RULE:FAIL") {
                    assert!(mtype != "SpliceSite");
                    fifty_bp_rule_fails_but_conserved += 1;
                }
                if mtype != "SpliceSite" {
                    assert!(info.starts_with("GERP_DIST:"), "{}", line);
                } else {
                    assert!(info.starts_with("INTRON_SIZE:"), "{}", line);
                }
            }
            let alt = alt.replace('N', "A");
            observed_input.push_str(&format!("{} {} {} {}\n", chr, pos, reference, alt));
            expected.insert(
                (pos.to_string(), alt),
                (mtype.to_string(), label.to_string(), info.to_string()),
            );
        }
        let observed_file = dir.join(format!("{}.observed.txt", id));
        fs::write(&observed_file, observed_input).unwrap();
        let classified_file = dir.join(format!("{}.classified.tsv", id));
        let mut args = vec![
            "--action",
            "classify",
            "--gff3",
            GFF3,
            "--genome",
            GENOME,
            "--observed-mutations",
            p(&observed_file),
            "--classified-mutations",
            p(&classified_file),
            "--id",
            id,
        ];
        args.extend(loftee);
        stdout(&args);
        let classified = fs::read_to_string(&classified_file).unwrap();
        let mut lines = classified.lines();
        assert_eq!(
            lines.next().unwrap(),
            "region\tchromosome\tposition\tmutation_type\tchange\tloftee\tloftee_info"
        );
        let mut differences = 0;
        let mut compared = 0;
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            let pos = (f[2].parse::<usize>().unwrap() + 1).to_string();
            let alt = f[4].split("->").nth(1).unwrap().to_string();
            let observed = (type_name(f[3]).to_string(), f[5].to_string(), f[6].to_string());
            let e = &expected[&(pos, alt)];
            // classify and enumerate may differ in type for a few indels at
            // the CDS start (an existing difference); compare the rest
            if e.0 == observed.0 {
                compared += 1;
                if *e != observed {
                    differences += 1;
                    eprintln!("{}: possible_mutations {:?}, classify {:?}", id, e, observed);
                }
            }
        }
        assert!(compared > 1000, "{}: only {} compared", id, compared);
        assert_eq!(differences, 0, "{}", id);
    }
    eprintln!("{:?}", labels);
    for label in [
        "Nonsense HC",
        "Nonsense LC:END_TRUNC",
        "FrameshiftIndel HC",
        "FrameshiftIndel LC:END_TRUNC",
        "SpliceSite HC",
    ] {
        assert!(labels.contains_key(label), "no {} in {:?}", label, labels);
    }
    assert!(fifty_bp_rule_fails_but_conserved > 0);
}

/// Sum of the probabilities of `possible_mutations` (column 8) per type and
/// LOFTEE label (column 11, if any).
fn sums(sites: &str, with_labels: bool) -> BTreeMap<(String, String), f64> {
    let mut out = BTreeMap::new();
    for line in sites.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f.len(), if with_labels { 12 } else { 10 }, "{}", line);
        let label = if with_labels {
            f[10].split(':').next().unwrap().to_string()
        } else {
            ".".to_string()
        };
        *out.entry((f[5].to_string(), label)).or_insert(0.0) += f[7].parse::<f64>().unwrap();
    }
    out
}

fn expected_counts(dir: &Path, extra: &[&str]) -> BTreeMap<String, f64> {
    let (pos, exp) = (dir.join("possible.txt"), dir.join("expected.tsv"));
    let mut args = vec![
        "--action",
        "enumerate",
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--id",
        IDS[0],
        "--possible-mutations",
        p(&pos),
    ];
    args.extend(model());
    args.extend(extra);
    stdout(&args);
    stdout(&[
        "--action",
        "expect",
        "--possible-mutations",
        p(&pos),
        "--expected-mutations",
        p(&exp),
        "--id",
        IDS[0],
    ]);
    let text = fs::read_to_string(exp).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let values: Vec<&str> = lines.next().unwrap().split('\t').collect();
    header
        .iter()
        .zip(&values)
        .skip(1)
        .map(|(h, v)| (h.to_string(), v.parse().unwrap()))
        .collect()
}

fn model() -> Vec<&'static str> {
    vec![
        "--point-mutation-probabilities",
        "SNV_mutation_rate_model.txt",
        "--indel-mutation-probabilities",
        "indel_mutation_rate_model.txt",
        "--scaling-factor",
        "8586",
    ]
}

/// `--loftee` only adds labels; `--loftee-filter` leaves out the LC mutations
/// in enumerate, possible_mutations and classify.
#[test]
fn filter_counts_high_confidence_lof_only() {
    let dir = out_dir("genovo_loftee_filter");
    let gerp = write_gerp(&dir);
    let site_args = |extra: &[&str]| {
        let mut a = vec![
            "--action",
            "possible_mutations",
            "--gff3",
            GFF3,
            "--genome",
            GENOME,
            "--id",
            IDS[0],
        ];
        a.extend(model());
        a.extend(extra);
        a.iter().map(|s| s.to_string()).collect::<Vec<String>>()
    };
    let as_str = |v: &Vec<String>| v.iter().map(String::as_str).collect::<Vec<&str>>().join(" ");
    let plain = stdout(&site_args(&[]).iter().map(String::as_str).collect::<Vec<_>>());
    let labelled_args = site_args(&["--loftee", "--gerp", p(&gerp)]);
    let labelled = stdout(&labelled_args.iter().map(String::as_str).collect::<Vec<_>>());
    // the same mutations, types and probabilities
    let first10 = |text: &str| -> Vec<String> {
        text.lines()
            .map(|l| l.split('\t').take(10).collect::<Vec<_>>().join("\t"))
            .collect()
    };
    assert_eq!(first10(&plain), first10(&labelled), "{}", as_str(&labelled_args));
    let by_label = sums(&labelled, true);
    let hc = |t: &str| {
        by_label
            .get(&(t.to_string(), "HC".to_string()))
            .copied()
            .unwrap_or(0.0)
    };
    let lc = |t: &str| {
        by_label
            .get(&(t.to_string(), "LC".to_string()))
            .copied()
            .unwrap_or(0.0)
    };
    assert!(lc("Nonsense") > 0.0 && hc("Nonsense") > 0.0);
    assert!(lc("FrameshiftIndel") > 0.0 && hc("FrameshiftIndel") > 0.0);

    // enumerate + expect: LoF counts are the HC sums
    let filter = ["--loftee-filter", "--gerp", p(&gerp)];
    let all = expected_counts(&out_dir("genovo_loftee_filter_all"), &[]);
    let filtered = expected_counts(&out_dir("genovo_loftee_filter_hc"), &filter);
    let close = |a: f64, b: f64| (a - b).abs() <= 1e-5 * a.abs().max(1e-12);
    for (name, t) in [
        ("nonsense", "Nonsense"),
        ("frameshift_indel", "FrameshiftIndel"),
        ("splice_site", "SpliceSite"),
    ] {
        assert!(
            close(filtered[name], hc(t)),
            "{}: {} vs HC {}",
            name,
            filtered[name],
            hc(t)
        );
        assert!(close(all[name], hc(t) + lc(t)), "{}", name);
    }
    for name in ["synonymous", "missense", "inframe_indel"] {
        assert!(close(filtered[name], all[name]), "{}", name);
    }
    // possible_mutations with filtering: no LC lines, the same sums
    let filtered_sites = stdout(&site_args(&filter).iter().map(String::as_str).collect::<Vec<_>>());
    assert!(filtered_sites
        .lines()
        .all(|l| !l.split('\t').nth(10).unwrap().starts_with("LC")));
    assert!(close(
        sums(&filtered_sites, true)[&("Nonsense".to_string(), "HC".to_string())],
        hc("Nonsense")
    ));
    // with --include-unknown, the LC mutations are written as Unknown
    let with_unknown_args = site_args(&["--loftee-filter", "--gerp", p(&gerp), "--include-unknown"]);
    let with_unknown = stdout(&with_unknown_args.iter().map(String::as_str).collect::<Vec<_>>());
    let lc_lines: Vec<&str> = with_unknown
        .lines()
        .filter(|l| l.split('\t').nth(10).unwrap().starts_with("LC"))
        .collect();
    assert!(!lc_lines.is_empty() && lc_lines.iter().all(|l| l.split('\t').nth(5) == Some("Unknown")));

    // classify: LC mutations are Unknown (code 0)
    let cls = dir.join("classified.tsv");
    let mut args = vec![
        "--action",
        "classify",
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--observed-mutations",
        OBSERVED,
        "--classified-mutations",
        p(&cls),
    ];
    args.extend(filter);
    stdout(&args);
    let classified = fs::read_to_string(&cls).unwrap();
    let mut n_lc = 0;
    for line in classified.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        if f[5].starts_with("LC") {
            assert_eq!(f[3], "0", "{}", line);
            n_lc += 1;
        } else if f[5] == "HC" {
            assert!(matches!(f[3], "3" | "6" | "9"), "{}", line);
        }
    }
    eprintln!("{} observed LC mutations", n_lc);
}

#[test]
fn option_errors() {
    let o = run(&[
        "--action",
        "possible_mutations",
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--id",
        IDS[0],
        "--loftee",
        "--filter-plof",
    ]);
    assert!(!o.status.success());
    let o = run(&[
        "--action",
        "possible_mutations",
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--id",
        IDS[0],
        "--loftee",
        "--loftee-gerp-cutoff",
        "x",
    ]);
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("--loftee-gerp-cutoff"));
    // negative cutoffs are numbers, without --gerp a warning
    let o = run(&[
        "--action",
        "possible_mutations",
        "--gff3",
        GFF3,
        "--genome",
        GENOME,
        "--id",
        IDS[0],
        "--loftee",
        "--loftee-gerp-cutoff",
        "-58",
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(String::from_utf8_lossy(&o.stderr).contains("END_TRUNC"));
    let text = String::from_utf8(o.stdout).unwrap();
    assert!(text.lines().any(|l| l.contains("\tHC\tBP_DIST:")));
    assert!(!text.contains("GERP_DIST"));
}
