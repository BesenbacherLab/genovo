# GeNovo - Identifying disease genes with de-novo mutations

Everyone's genome contains some de-novo mutations: Mutations that are not part
of their parents' genomes. These de-novo mutations can cause diseases but it is
still difficult to tell the disease-causing mutations apart from the ones that
have nothing to do with the disease. Furthermore, some genes have a higher
chance of mutating which means that a mutation in that gene is not that
surprising, even when that gene is mutated in several people with the disease.

The name GeNovo is a wordplay of *gene* detection based on *de-novo* mutations.

## Installation

First make sure you have `cargo` installed (`cargo` is the package manager for the Rust programming language).
The easiest way is to install [rustup](https://rustup.rs/) which automatically installs `cargo` and `rustc` for you.
If you are on a UNIX-based system, you might also be able to install `cargo` or `rustup` with your system's package manager.

Secondly, you need to run

```
cargo install genovo
```

This will install the program genovo at a specific location which will be shown at the end of the installation process.
(You might need to make adjustments to your `PATH` variable afterwards to automatically find the program on the command line)

To build from this repository (with `mutexpect` and `pattern_partition_prediction` checked out next to it), use the dependency versions in `Cargo.lock`, which the tests were run with:

```
cargo build --release --locked
```

This needs Rust 1.75 or newer (`rustc --version`; update with `rustup update stable`). The binary is `target/release/genovo`.

The pipeline is a single executable (`genovo`) that combines the important steps from an earlier pipeline that I had put together.
(A single executable is much easier for others to use than a bunch of scripts that are all loosely coupled.)

## The pipeline

![The pipeline flow diagram](https://raw.githubusercontent.com/jbethune/kmerPaPa/master/images/pipeline.png)

* Input files are marked with a blue border.
* Pipeline steps are marked with a black border and their names are marked with a green background
* Optional parameters are marked in blue font

### Input and output files

The pipeline takes several parameters that represent input and output parameters.
Depending on which action is being performed, a parameter is treated either as an input or output file.

For example, the `--possible-mutations` parameter is treated as an input file when `--action expect` is specified.
But it is treated as an output file when `--action enumerate` is specified.
If no `--action` parameter is specified, the entire pipeline is run and all file parameters are treated as output files.

You can also specify `-` as an input/output file. Depending on the action, this will be interpreted as STDIN or as STDOUT.

Steps that need the genomic regions accept either `--genomic-regions` (the output of `transform`) or `--gff3` (then the transformation is done on the fly).

Progress messages (the name of each step as it starts) are written to STDERR, so STDOUT only contains results.

If you are concerned about space usage you can also specify a `.gz` file extension.
`genovo` will automatically decompress/compress the input/output.

### ID filtering / parallel processing

If the parameter `--id ID` is passed, then only the transcripts with the id `ID` will be processed
and all other transcripts will be ignored.
This is useful for parallel execution. You can launch a separate job for each transcript.
Please note that `--id` only accepts a single ID. Multiple IDs are not supported.
You can always combine the `--id` option with any other option.


### Running the whole pipeline on a test data set

You can download a test data set from the git repository:
```
wget https://github.com/BesenbacherLab/genovo/raw/master/test_data/test_data.tar.gz
tar xvzf test_data.tar.gz
```

If you want to run the entire pipeline on the test data you can execute the following command:

```
genovo \
	--gff3 gencode.v19.CCDS.annotation.protein_coding_chr22.gff3.gz \
	--observed-mutations observed_mutations_DDD_2017_chr22.txt \
	--genome hg19_chr22.2bit \
	--point-mutation-probabilities SNV_mutation_rate_model.txt \
	--indel-mutation-probabilities indel_mutation_rate_model.txt \
	--scaling-factor 8586 \
	--significant-mutations genovo_results_chr22.txt
```

The command takes ~1 minute to run on a "normal" computer.
It will write the comparisons between observed, expected and sampled mutations to the file `genovo_results_chr22.txt`. If the `--significant-mutations` argument is not used the results will be written to STDOUT.
For each transcript, *t*, and mutation type, *m*, the output contains the following output columns

|column           |description|
|---|---|
|observed        |The observed number of *m* mutations in transcript *t*|
|expected        |The expected number of *m* mutations in transcript *t* according to the mutation rate models|
|expected_lower  |The lower bound of the confidence interval around the expected number of mutations|
|expected_upper  |The upper bound of the confidence interval around the expected number of mutations|
|p_value         |A one-sided p-value for a sampling based test of whether there are more observed mutations of type *m* in transcript *t* than we would expect given our mutation rate models|


When testing on de novo mutations from trios the `--scaling-factor` should be two times the number of trios in the data. So in this case it is set to 8586 since the data set with the observed mutations looked for de novo mutations in 4293 children.

The `SNV_mutation_rate_model.txt` and `indel_mutation_rate_model.txt` are mutation rate models created using [kmerpapa](https://github.com/BesenbacherLab/kmerPaPa). Trained models in the right input format can be downloaded from https://github.com/BesenbacherLab/Genovo_Input.


### Steps of the pipeline

Each of the steps of the pipeline can be executed in isolation by specifying the `--action STEP` parameter.
If no `--action` parameter is specified, all steps are run (as in the example above).

#### transform

This step takes a [gff3 file](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#readme) and transforms it into a `regions` file (which is a file format that is specific to genovo).
The gff3 file needs to have the following attributes:

* ID: For identifying the name of each transcript
* Parent: For identifying the parent of exons and CDS regions

The file must also be sorted so that all exons and CDS regions are listed after their transcript line.

Some example gff3 files can be obtained from [GENCODE](https://www.gencodegenes.org/human/)

It is highly recommended that you only use annotations that contain exclusively coding transcripts.
The results for non-coding transcripts are a lot less informative.

##### parameters

|parameter        |input/output|description|
|---|---|---|
|--gff3           |input |The GFF3 annotation file that contains transcript, exon and CDS information             |
|--genomic-regions|output|A transformed version of the GFF3 entries, limited to transcripts, exons and CDS regions|
|--regions-bed    |output (optional)|The positions where per-site allele numbers are needed (exons ± 2 bases; whole transcripts with `--include-intronic`) as a BED file. Can also be made from an existing `--genomic-regions` file.|


#### enumerate

This goes through the genomic sequence of every transcript and determines all possible point mutations together with their probability
This command takes a [PaPa-rates file](https://crates.io/crates/pattern_partition_prediction) which represents the point mutation probabilities for different sequence contexts as described in our paper *Improved prediction of site-specific mutation rates using k-mer pattern partition *.
This command also takes a reference genome file in the [twobit file format](https://genome.ucsc.edu/goldenPath/help/twoBit.html). Make sure that all your input files use the same genome build!

The output file is fairly simple:
```
#transcript_ID
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
#other_transcript_ID
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
numeric_mutation_type:mutation_probability
```

The possible mutations for each transcript start with a `#ID` line followed
with one line for every possible mutation where the probability is not `NaN`.
Keep in mind that every position may mutate into one of the 3 other bases, so
that you have 3 possible point mutations for every position.

The numeric mutation type is defined in the section [Mutation type numeric code table](#mutation-type-numeric-code-table).
 
##### parameters

|parameter                     |input/output|description|
|---|---|---|
|--genomic-regions             |input  |A transformed version of the GFF3 entries, limited to transcripts, exons and CDS regions|
|--genome                      |input  |A twobit reference genome sequence file|
|--point-mutation-probabilities|input  |A pattern partition file with the sequence-context-dependent point mutation probabilities|
|--indel-mutation-probabilities|input  |A pattern partition file with the sequence-context-dependent indel probabilities|
|--possible-mutations          |output |The possible mutations for each transcript|
|--scaling-factor              |number |Scaling of the per-generation probabilities r to 1 − (1 − r)^s (default 1), e.g. 2 × the number of trios for de novo mutations|
|--allele-numbers              |input (optional)|Per-site allele numbers, see [Uneven coverage](#uneven-coverage-per-site-allele-numbers)|
|--calibration                 |input (optional)|Instead of `--scaling-factor`: scaling that depends on the allele number|
|--min-an                      |number (optional)|Leave out sites with a smaller allele number|
|--loftee-filter, --gerp       |input (optional)|Count only LOFTEE's high-confidence LoF mutations, see [LoF confidence](#lof-confidence-loftees-filters)|


#### classify

This step classifies observed point mutations (which is a raw input file you should provide).
The input file format is as follows:

```
chr1 1230448 G A
chr1 1609723 C T
chr1 1903276 C T
chr1 2574999 C T
```

The first column is the chromosome (must be the same naming convention as in the annotations and the twobit file).
The second column is the position on the chromosome. This position has to be 1-based.
The third column is the reference base. If this reference base does not match what is in the twobit file, you will get an error.
The fourth column is the new base that was mutated into.

##### parameters

|parameter             |input/output|description|
|---|---|---|
|--observed-mutations  |input |Observed point mutations from your own data|
|--genome              |input |A twobit reference genome sequence file|
|--genomic-regions     |input |A transformed version of the GFF3 entries, limited to transcripts, exons and CDS regions|
|--classified-mutations|output|The observed mutatations. But annotated with numeric mutation types|
|--allele-numbers, --min-an|input (optional)|Leave out observed mutations at the sites that `enumerate` leaves out (give the same values to both steps)|
|--loftee or --loftee-filter, --gerp|input (optional)|LOFTEE's labels (columns `loftee`, `loftee_info`), see [LoF confidence](#lof-confidence-loftees-filters)|

#### expect

This step sums up all possible mutations and their probabilities and gives you
a grand total that represents the expected number of mutations for each gene.

##### parameters

|parameter            |input/output|description|
|---|---|---|
|--possible-mutations |input |The possible mutations for each transcript|
|--expected-mutations |output|The number of expected mutations per transcript for each mutation type|

#### sample

This step randomly picks possible mutations and simulates randomly if they occur based on their probability.
This creates an empirical distribution of sampled mutations that will be used to calculate a p-value in a later step.

##### parameters

|parameter                 |input/output|description|
|---|---|---|
|--possible-mutations      |input  |The possible mutations for each transcript|
|--number-of-random-samples|integer|How many random samples? default=1000|
|--sampled-mutations       |output |Distributions of sampled mutations|
|--seed                    |integer|Seed for the random numbers (optional, see below)|

##### Reproducibility

The `sample` and `compare` steps use random numbers. With `--seed N` a run is reproducible: the same inputs and seed give identical output files.
Each transcript and mutation type gets its own random number stream derived from the seed, so the results are the same whether the pipeline is run in one go or step by step, and whether all transcripts are processed together or one at a time with `--id`.
Without `--seed`, a random seed is chosen and written to STDERR, so the run can still be repeated.
Give the same seed to `sample` and `compare` when running the steps separately.

#### compare

Compare the observed, expected and sampled mutations to determine a p-value for
observing that many mutations in a transcript by chance.
Please note that the output file contains the results of all statistical tests, even the non-significant ones.
And note that the p-values does not take into account the expected number of mutations.

##### parameters

|parameter               |input/output|description|
|---|---|---|
|--classified-mutations  |input |The observed mutatations. But annotated with numeric mutation types|
|--expected-mutations    |input |The number of expected mutations per transcript for each mutation type|
|--sampled-mutations     |input |Distributions of sampled mutations|
|--significant-mutations |output|A table of observed, expected, sampled and a p-value for different mutation types.|
|--seed                  |integer|Seed for the random pairing of samples when combining LoF types (optional, see `sample`)|


#### Uneven coverage: per-site allele numbers

For population data with uneven coverage (e.g. exomes), the probability of
observing a mutation depends on how many chromosomes were called at its site.
`enumerate` then computes the probability of a mutation with per-generation
rate r at a site with allele number n as 1 − (1 − r)^(S · h(n)), and
`classify` leaves out the observed mutations at the same masked sites. The
design is described in `docs/design/allele-number-scaling.md`.

* `--allele-numbers FILE`: the number of called chromosomes at each site,
  tab-separated, 1-based, sorted by chromosome and position, either as
  `chrom pos an` or as `chrom:pos an` (the layout of the gnomAD v4.1
  all-sites AN files), plain or gzip/bgzip-compressed, with an optional header.
  Sites that are not in the file are not callable and are left out of both
  expected and observed counts. For indels, the allele number of the base
  after the insertion/deletion point is used. To combine several data sets
  (e.g. gnomAD exomes and genomes), add up their allele numbers per site
  beforehand.
* `--calibration FILE`: the scaling factors `S` (separately for point
  mutations and indels) and the relative effective sample size `h(n)`
  (linear between the listed allele numbers, constant outside them), in this
  format:

  ```
  # genovo calibration v1
  # S_snv	2.2e7
  # S_indel	1.9e7
  # n_min	50000
  an	h
  50000	0.081
  1612000	1
  ```

  Without `--calibration`, `S = --scaling-factor` and `h = 1`.
* `--min-an N`: also leave out sites with an allele number below `N`
  (overrides `n_min` of the calibration).

**Fitting the calibration** (`--action calibrate`): `S` and `h` are fitted on
the possible synonymous point mutations (each site and alternative base once)
and the observed variants (`--observed-mutations`, e.g. all variants of the
population data set; indels are ignored):

```
genovo --action calibrate --gff3 annotation.gff3.gz --genome genome.2bit \
    --point-mutation-probabilities SNV_mutation_rate_model.txt \
    --allele-numbers allele_numbers.tsv.gz --observed-mutations variants.txt \
    --calibration calibration.tsv [--min-an N] [--calibration-an-step 1.05] \
    [--calibration-min-observed 1000] [--calibration-monotone] \
    [--calibration-ids canonical_transcripts.txt]
```

The sites are split into bins of allele numbers: the allele numbers are cut
in steps of 5% (`--calibration-an-step`), and neighbouring steps are joined,
from the highest allele number down, until a bin has at least 1000 observed
and 1000 unobserved possible synonymous mutations (`--calibration-min-observed`).
Allele numbers that only few sites have (e.g. sites covered only by the genomes
of a data set that combines exomes and genomes) thus get bins of their own. In
each bin, the scaling is chosen so that the expected number of observed
variants equals the observed number; `S` is the scaling of the best-covered
bin, and `h` of each bin is its scaling divided by `S` (constant within the
bin, linear between bins). `h` may decrease with the allele number and exceed 1, because
the allele number also reflects the quality of the data at a site, not only
the number of samples. With `--calibration-monotone`, bins are instead merged
until `h` does not decrease (and `h` is at most 1).
`S_indel` is not fitted (it is set to `S_snv`; edit the file if needed).
Diagnostics (observed and expected counts with the fitted calibration, by
allele-number bin and by mutation-rate decile, for all sites and separately
for CpG transitions, other transitions and transversions) are written to
`<calibration>.diagnostics.tsv` (or `--calibration-diagnostics FILE`); the
observed/expected ratio should be close to 1 everywhere.

The allele number file is read once per run, keeping only the sites in the
transcripts. A genome-wide file takes several minutes to read, so for many runs
(e.g. one per transcript with `--id`) restrict it beforehand to the BED file
written by `--action transform --regions-bed FILE`.

#### LoF confidence: LOFTEE's filters

gnomAD counts only the putative loss-of-function (LoF) variants that
[LOFTEE](https://github.com/konradjk/loftee) rates high confidence (HC).
genovo applies the LOFTEE filters that gnomAD v4.1 uses (without LOFTEE's
splice predictions) to nonsense, frameshift and splice-site mutations
(`src/loftee.rs`):

* `END_TRUNC` (nonsense, frameshift): at most 50 bases upstream of the last
  exon junction as LOFTEE measures it, and the GERP++ scores from the mutation
  to the stop codon sum to at most `--loftee-gerp-cutoff` (0 as in gnomAD
  v4.1.1; LOFTEE's default is -58). Needs `--gerp`. As in LOFTEE's code, the
  exons downstream of the stop codon (3' UTR exons) count towards that
  distance, so `END_TRUNC` rarely applies in transcripts with a 3' UTR intron.
* `SMALL_INTRON`, `5UTR_SPLICE`, `3UTR_SPLICE` (splice site): an intron
  shorter than 15 bases, or upstream of the start codon or downstream of the
  stop codon.

Options:

* `--loftee`: label the LoF mutations (`HC`, or `LC:` and the filters) in the
  output of `possible_mutations` (two more columns) and `classify` (columns
  `loftee` and `loftee_info`), with the numbers the filters use; the mutation
  types are unchanged.
* `--loftee-filter`: as `--loftee`, and count the low-confidence (LC) mutations
  as `Unknown` in `enumerate`, `possible_mutations` and `classify`, i.e. LoF
  counts are high confidence only. Conflicts with `--filter-plof`.
* `--gerp FILE`: GERP++ rejected-substitution scores per site (layout as
  `--allele-numbers`: `chrom pos score` or `chrom:pos score`, 1-based, sorted;
  only the exon sites are read), e.g. from LOFTEE's
  `gerp_conservation_scores.homo_sapiens.GRCh38.bw` at the sites of
  `--action transform --regions-bed`. LOFTEE takes the mean scores of the
  intervals from the bigWig summary levels, which approximate them, so sums
  close to the cutoff may come out differently.

Not applied: `GC_TO_GT_DONOR` (genovo's splice sites have a GT donor),
`ANC_ALLELE` (needs the human ancestral sequence) and LOFTEE's flags.

#### possible_mutations

Not part of the pipeline (it only runs when given with `--action possible_mutations`).
Writes every possible mutation of every transcript to STDOUT, one tab-separated line per mutation:

|column|description|
|---|---|
|1|chromosome|
|2|position (1-based); for indels the anchor base before the insertion/deletion point, as in VCF and in `--observed-mutations`|
|3|reference base (for indels: the anchor base)|
|4|alternative: a base for a point mutation; the anchor base followed by `NNN` for an in-frame indel or by `N` for a frameshift indel|
|5|transcript ID|
|6|mutation type (e.g. `Missense`, `SpliceSite`, `FrameshiftIndel`)|
|7|sequence context centred on the mutated base (for indels: on the base after the insertion/deletion point)|
|8|only with `--point-mutation-probabilities`: the probability of the mutation, as used by `enumerate`|
|9|only with `--point-mutation-probabilities`: the allele number of the site (`NA` without `--allele-numbers`)|
|10|only with `--point-mutation-probabilities`: the per-generation rate of the model, before scaling|
|last two|only with `--loftee` or `--loftee-filter`: LOFTEE's label of LoF mutations (`HC`, or `LC:` and the filters; `.` for other mutations) and the numbers the filters are based on (`GERP_DIST`, `BP_DIST`, `DIST_FROM_LAST_EXON`, `50_BP_RULE`, `INTRON_SIZE`, as in LOFTEE's `LoF_info`), see [LoF confidence](#lof-confidence-loftees-filters)|

The same mutations are considered as in `enumerate`, with the same filters (`--include-intronic`, `--include-unknown`, `--filter-plof`, `--id`).
The sequence context has the k-mer size of the `--point-mutation-probabilities` model if one is given, and is 9 bases otherwise.

With `--point-mutation-probabilities` (and `--indel-mutation-probabilities` for indels), columns 8–10 are written; column 8 holds the probability of each mutation after `--scaling-factor` or `--calibration` and `--allele-numbers`, and the mutations that `enumerate` leaves out are not written (masked sites, undefined rates, indels without an indel model). Summed per transcript and mutation type, the probabilities are the expected numbers of `expect`, so the expected number of mutations in any subset of sites (e.g. by position in the coding sequence) is the sum over those sites.

```
genovo --action possible_mutations --gff3 annotation.gff3.gz --genome hg19.2bit > possible_mutations.tsv
```

##### parameters

|parameter        |input/output|description|
|---|---|---|
|--genomic-regions or --gff3|input|The transcripts|
|--genome         |input |The reference genome in 2bit format|
|--point-mutation-probabilities|input (optional)|The size of the sequence context, and columns 8–10|
|--indel-mutation-probabilities|input (optional)|Columns 8–10 for indels|
|--scaling-factor, --calibration, --allele-numbers, --min-an|input (optional)|Scaling and masking of the probabilities, as in `enumerate`|
|--loftee or --loftee-filter, --gerp, --loftee-gerp-cutoff|input (optional)|LOFTEE's labels, see [LoF confidence](#lof-confidence-loftees-filters)|


## Tests

`cargo test` runs the unit tests and end-to-end tests on three transcripts of the chr22 test data (`test_data/test_data.tar.gz`, unpacked automatically; needs `tar`).
The end-to-end tests check that observed mutations are classified like the corresponding possible mutations, and compare the pipeline output (with `--seed 1`, in one go and step by step) with stored results in `tests/data/`.
The whole of chr22 is tested with `cargo test --release -- --ignored`.
After an intended change of the results, update the stored files with `GENOVO_BLESS=1 cargo test` (and `GENOVO_BLESS=1 cargo test --release -- --ignored pipeline_chr22`) and review the diff.

## Mutation type numeric code table

The input and output files of the pipeline are almost all text files. In order
to be more space-efficient, there is a numeric encoding for the different types
of mutations.

| mutation type name            | numeric code |
|-------------------------------|--------------|
| unknown                       | 0            |
| synonymous                    | 1            |
| missense                      | 2            |
| nonsense                      | 3            |
| loss of stop codon            | 4            |
| loss of start codon           | 5            |
| loss of canonical splice site | 6            |
| Mutation within an intron     | 7            |
| Non-frameshift indel          | 8            |
| Frameshift indel              | 9            |
|                               |              |
