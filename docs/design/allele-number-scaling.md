# Design note: per-site allele numbers and scaling-factor calibration

Status: proposal, 2026-10-11. Implemented so far (§8): steps 1–3
(2026-10-16): `src/allele_numbers.rs`, `--regions-bed`, `src/scaling.rs`,
and `--allele-numbers`, `--calibration`, `--min-an` for `enumerate` and
`classify` (tests in `tests/allele_numbers.rs`).

## 1. Goal

Make genovo's expected counts, sampled distributions, p-values and confidence
intervals correct for population data with uneven coverage (gnomAD exomes, or
exomes and genomes combined). Today this needs a correction outside genovo.
The scaling should:

- use the number of called chromosomes **at each site** (the allele number,
  AN), not one number for the whole data set or one per transcript;
- be **fitted by genovo** from the data (`--action calibrate`) instead of
  tuned by hand;
- be applied **before** the saturation formula, so that well-saturated sites
  (CpG transitions) and low-rate sites respond correctly to low coverage;
- reduce to today's behaviour when no AN is given (de novo data, uniform WGS);
- not depend on a particular data source: genovo reads **one per-site count
  file** in a simple format. Anything source-specific, such as adding up the
  gnomAD exome and genome allele numbers, is done in a preprocessing step
  outside genovo (§4.1).

## 2. Current method and its problems

Today, `enumerate` turns each per-generation rate `r` into

```text
p = 1 − (1 − r)^S          S = --scaling-factor (2 × number of generations)
```

and `sample` draws every possible mutation once with probability `p`. For de
novo data `S = 2 × trios`. For gnomAD in the paper (Methods, "Scaling mutation
rates" and "Correction for coverage"):

1. `S` was chosen so that the expected synonymous count matched the observed
   count in the whole exome.
2. Each transcript's expected counts were divided by a LOWESS fit of
   synonymous O/E against the transcript's mean AN (Supplementary Fig. 6),
   where the mean AN came from the gnomAD v2 VCF, i.e. from variant sites only.
3. Transcripts with mean AN < 50,000 were removed.

Problems:

- **Applied after saturation.** Halving the number of called chromosomes
  roughly halves `p` at a low-rate site, but hardly changes a CpG transition
  with `p` close to 1. A single factor per transcript treats both alike, so
  CpG-rich categories (LoF via CGA→TGA) are over-corrected relative to
  synonymous.
- **Coverage varies within a transcript.** Exon ends, where splice sites are,
  are covered worse than exon middles. A transcript mean fitted on synonymous
  sites does not apply to the splice-site or LoF sites of the same transcript.
- **AN from variant sites only.** Variants are called more easily where
  coverage is good, so the mean AN is biased upwards, most for poorly
  covered transcripts.
- **Only the expectation is corrected.** `sample` and `compare` use the
  uncorrected probabilities, so p-values (in particular for depletion) and the
  sampled confidence intervals do not match the corrected expectations.
- **Whole transcripts are dropped** at AN < 50,000, even if most of their
  sites are well covered.

## 3. Model

For possible mutation `i` (a site and an alternative allele) with
per-generation rate `r_i` at a site with allele number `n_i`:

```text
p_i = 1 − (1 − r_i)^(S · h(n_i))           computed as −expm1(S · h(n_i) · ln1p(−r_i))
```

- `S`: the scaling factor of a fully covered site; separate values `S_snv`
  for point mutations and `S_indel` for indels (§6.3).
- `h(n)`: the relative effective sample size at a site with `n` called
  chromosomes, monotone non-decreasing, `h(n_ref) = 1` at a reference AN
  `n_ref` (default: the 99th percentile of `n` over the calibration sites).
  `h` is estimated from the data (§5), not assumed to be `n / n_ref`: the
  number of segregating sites grows sublinearly with sample size, and the
  mixture of capture kits gives plateaus such as the one at 150–180k in
  Supplementary Fig. 6.
- For an indel, `n_i` is the allele number of the base after the
  insertion/deletion point (the position at which `enumerate` lists it; the
  VCF anchor + 1), so that `enumerate` and `classify` use the same site.
- `n_i = 0` (site not callable) gives `p_i = 0`. Sites with `n_i < n_min`
  (`--min-an`) are **masked**: they contribute neither expected nor observed
  mutations (§4.2).

Without an AN file, `h ≡ 1` and `S_snv = S_indel = --scaling-factor`, which is
exactly today's behaviour. The sampling step is unchanged: each possible
mutation is drawn once with probability `p_i`, which matches gnomAD, where a
variant site counts once however many carriers it has.

### 3.1 Combining data sets

genovo works with a single allele number per site. If the observed variants
come from several data sets, the counts are combined in preprocessing. For the
gnomAD v4.1 joint (exome + genome) variant set, that means summing the exome
and genome AN at each site. Summing assumes that one called chromosome counts
the same in both data sets; whether that holds is checked with the
diagnostics (§7) and not modelled inside genovo. If exomes and genomes turn
out to differ, they can be analysed separately (each with its own AN file and
calibration), or the preprocessing can weight one of them.

## 4. Data and inputs

### 4.1 Per-site allele number file

genovo reads one file with the number of called chromosomes at each site
(`--allele-numbers FILE`):

```text
chrom   pos     an
chr1    10001   16
chr1    10002   78
```

- Tab-separated, 1-based positions, sorted by chromosome and position. The
  locus may also be written as one column, `chr1:10001<TAB>16` (the layout of
  the gnomAD all-sites AN files). An optional header line and `#` comment
  lines are skipped. Plain text or gzip/bgzip (`flate2::MultiGzDecoder`).
- Unsorted input (a position not after the previous one, or a chromosome
  appearing again later) is an error.
- Sites missing from the file have `an = 0`: not callable, and always masked
  when an allele number file is given (also with `n_min = 0`).
- Chromosome names must match the 2bit genome and the regions file.
- genovo streams the file once and keeps only the positions where it
  enumerates possible mutations: exons, plus the splice-site bases in the
  introns, or the whole region with `--include-intronic`. These are stored
  per chromosome as a sorted vector of `(position, an)`, looked up by binary
  search. Tens of millions of sites for all protein-coding transcripts fit
  easily in memory.
- A genome-wide file (about 3 × 10⁹ lines for WGS) works, but reading it
  takes minutes every time (about 9 million lines per second, so 6–10
  minutes). When genovo is run as many jobs (e.g. one per
  transcript with `--id`), the file should be restricted to the regions in
  preprocessing. To make that easy, `--action transform` can also write the
  needed positions as a BED file (`--regions-bed FILE`: exons ± a margin), for
  use with `tabix -R` or `bedtools intersect`.

**Preprocessing for gnomAD v4.1** (outside genovo). The all-sites AN files
have the form `locus<TAB>AN` with `locus = chr1:10001`, one file each for
exomes and genomes. A preprocessing script:

1. restricts both files to the BED from `--regions-bed` (optional, for
   speed);
2. for the joint data set, sums the two AN values per position (a sorted
   merge of the two files; a position missing from one file counts as 0).

The `chr1:10001` layout can be read directly. For exomes only, or genomes
only, the gnomAD file can be given to genovo as it is.

gnomAD v4.1 is on GRCh38, so the regions, the 2bit genome and the observed
variants must be GRCh38 too (the paper used hg19/GENCODE v19 with gnomAD v2).

The same format can be produced from coverage files (e.g. the number of
samples with depth ≥ 20 at each site, × 2) for data sets without all-sites
AN.

### 4.2 Observed variants

The observed variants are given in genovo's observed-mutations format
(`chr pos ref alt`), one line per alternative allele (each alternative allele
is a separate possible mutation). They must come from the same samples as the
allele numbers: for gnomAD, PASS variants with AC > 0 in the data sets whose
AN was summed. `classify` gets `--allele-numbers` and `--min-an`, so
that observed variants at masked sites are dropped exactly where the expected
side drops them. The site's AN is always taken from the allele number file,
not from the variant's VCF record, so that both sides use the same number.

## 5. Calibration: `--action calibrate`

```text
genovo --action calibrate \
    --genomic-regions regions.txt --genome hg38.2bit \
    --point-mutation-probabilities snv_model.txt \
    --indel-mutation-probabilities indel_model.txt \
    --allele-numbers allele_numbers.tsv.gz \
    --observed-mutations gnomad_variants.txt \
    --calibration calibration.tsv            # output
```

The calibration needs each possible mutation's position (to look up `n_i` and
whether it was observed), which the possible-mutations file does not keep. So
`calibrate` walks the transcripts itself with `for_each_possible_mutation`,
as `enumerate` does.

### 5.1 Calibration sites

- Point mutations: the possible **synonymous** mutations. This is the usual
  near-neutral choice (gnomAD constraint uses it, and so did the paper).
  Overlapping transcripts share sites, so each `(position, alt)` is used once,
  optionally only from canonical / MANE transcripts (`--calibration-ids FILE`).
- Indels: see §6.3.
- Masked sites (`n_i < n_min`) are left out.
- Each unit has `(r_i, n_i, y_i)`, with `y_i = 1` if that exact alternative
  allele is observed.

### 5.2 Estimating `S` and `h`

`h` is a monotone step or piecewise-linear function on bins of `n`:

1. Split the calibration units into `B` bins of `n` (default 25 quantile bins,
   with the top bin holding `n ≥ n_ref`).
2. For each bin `b`, find the scaling `s_b` that makes the expected count
   equal the observed count:
   `Σ_{i∈b} [1 − (1 − r_i)^{s_b}] = Σ_{i∈b} y_i`.
   The left side increases monotonically with `s_b`, so bisection finds it.
   With rare, independent units this is the maximum-likelihood estimate of a
   one-parameter-per-bin model.
3. `S = s_top` (the best-covered bin) and `h_b = s_b / S`. Then enforce
   monotonicity with pool-adjacent-violators, weighted by the expected count
   per bin, and set `h = 0` below `n_min`.
4. Between bin centres, `h(n)` is linear interpolation; above the top centre
   it is constant 1.

This is the site-level version of what the paper did with transcripts and
LOWESS. Fitting `S` on the best-covered sites, and scaling everything else
down from there, is the same idea, now applied before saturation and per
site.

Alternatives, if the bins turn out to be noisy:

- a parametric `h(n) = (n / n_ref)^β`, with a likelihood fit of `S` and `β`;
- a joint maximum-likelihood fit of `S` and all `h_b`, with the Bernoulli
  likelihood `Σ y log p + (1 − y) log(1 − p)`.

The bin-matching estimator is simple, deterministic and easy to inspect, so
it is the default.

### 5.3 Output

`calibration.tsv` holds the model, plus the provenance needed to reproduce it
(the numbers below are only illustrative):

```text
# genovo calibration v1
# S_snv      2.2e7
# S_indel    1.9e7
# n_ref      1612000
# n_min      50000
# source     synonymous, 18,412,345 units, 3,112,876 observed
an        h
50000     0.081
...
1612000   1.000
```

The diagnostics (§7) are written next to it.

## 6. Using the calibration

### 6.1 `enumerate`

New options: `--calibration FILE --allele-numbers FILE`.

- For each possible mutation, look up `n_i`, compute `h(n_i)`, and apply
  `p_i = −expm1(S · h(n_i) · ln1p(−r_i))` with `S_snv` or `S_indel`.
- Masked sites produce no event.
- `--scaling-factor` and `--calibration` are mutually exclusive.

`expect`, `sample` and `compare` are unchanged: they only see the `p_i`, so
expected counts, sampled distributions, p-values and confidence intervals are
all consistent.

### 6.2 `classify`

With `--allele-numbers` and `--min-an`, observed mutations at masked
sites are dropped (and counted in a log message), so that observed and
expected counts cover exactly the same sites.

### 6.3 Indels

Indels need their own `S_indel`:

- the indel model's rates were scaled separately (0.68 × 10⁻⁹ per base per
  generation);
- gnomAD calls indels less sensitively than SNVs;
- exome coverage affects indel calling differently.

There is, however, no synonymous-like neutral class for coding indels.
Options, roughly in order of preference:

1. **Non-coding indels near exons.** Enumerate possible indels in UTRs and
   in introns beyond the splice region (bases 10–50 from exon ends, still
   inside the exome capture), for calibration only, and fit `S_indel` on them
   with the same `h` as SNVs.
2. **Ratio from genomes.** In gnomAD genomes (uniform coverage), fit
   `S_indel / S_snv` on non-coding indels and SNVs genome-wide, and apply that
   ratio to the exome `S_snv`.
3. **Fallback:** fit `S_indel` on in-frame indels, accepting that they are
   under some selection (this biases frameshift O/E upwards).

Option 1 needs `for_each_possible_mutation` to optionally produce non-coding
indels. It does not today, because indels are only enumerated inside CDSs.

## 7. Validation

`calibrate` writes, and the workflow plots:

1. **Flatness by coverage.** Synonymous observed/expected per bin of `n`, for
   all sites, CpG transitions, other transitions and transversions
   separately. All should be ≈ 1 across `n`. The CpG split tests the
   interaction between saturation and coverage that the transcript-level
   correction gets wrong.
2. **Calibration by rate.** Observed fraction against mean `p` in bins of
   `r`, separately for high and low `n`. This checks the saturation shape
   `1 − (1 − r)^s`; if it fails at high `r`, revisit the functional form.
3. **Gene level.** Synonymous O/E per transcript against mean `n`, as
   Supplementary Fig. 6 but after correction: no trend. Missense and LoF O/E
   should no longer correlate with coverage among genes that are not
   constrained.
4. **Held-out chromosomes.** Calibrate on even chromosomes and check flatness
   on odd ones.
5. **Comparison with the paper.** LOEUF / LoF O/E against gnomAD's values,
   compared with the transcript-level LOWESS correction: the number of usable
   transcripts and their correlation with gnomAD.

## 8. Implementation plan

1. Reader for the allele number file (§4.1), restricted to the regions while
   streaming, and `--regions-bed` for `transform`. Unit tests with a small
   synthetic file (missing sites, header, gzip/bgzip input, unsorted input
   rejected).
2. Site-AN lookup structure plus the `h`/`S` calibration type (load and save
   `calibration.tsv`), and `p_i` computed with `expm1`/`ln1p`. Done: the
   precise formula is already used for `--scaling-factor`; it changed 196 of
   8,140 expected values on chr22 in the last float32 digit (≤ 2·10⁻⁷
   relative), no sampled distribution or p-value. Unit tests:
   `h ≡ 1` reproduces `--scaling-factor` exactly; monotone interpolation;
   masking.
3. `enumerate` and `classify` options (§6). Test: the chr22 golden results
   are unchanged without AN; with a synthetic AN file, masked sites disappear
   from both sides.
4. `calibrate` action with diagnostics (§5, §7). Test on simulated data:
   draw `y_i` from a known `S` and `h`, check that they are recovered.
5. Workflow: download and prepare the AN files, calibrate, run genovo per
   transcript on the cluster, make the diagnostic plots.
6. Indel calibration (§6.3), after looking at the SNV results.

## 9. Open questions

- Is summing exome and genome AN adequate for the gnomAD joint data, or
  should exomes and genomes be analysed separately (§3.1)?
- Should the calibration sites be restricted to canonical / MANE transcripts?
- How should indels be calibrated (§6.3)?
- What about chrX (hemizygous males; possibly different mutation rate
  scaling) and chrY?
- Is a single `h` enough, or does the relation between coverage and
  detection differ by variant type (e.g. CpG vs non-CpG, indels) or by
  ancestry composition? The diagnostics in §7 will show.
- `n_min` default: the paper used a mean AN of 50,000 per transcript; per
  site, a lower threshold may be fine, because the expectation is corrected
  rather than only filtered.
