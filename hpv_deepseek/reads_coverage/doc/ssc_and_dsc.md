# SSC and DSC UMI Consensus Pipeline on ERIS Nucleus

`ssc_and_dsc.py` generates and optionally submits SLURM jobs for single-strand consensus (SSC) and duplex consensus (DSC) UMI processing from paired-end FASTQ samples in one run folder.

The original shell workflow is preserved in `ssc_and_dsc.sh`. The Python version targets ERIS Nucleus and writes a master analysis log before submission.

## Usage

```bash
python ssc_and_dsc.py \
    --run-id <run_id> \
    [options]
```

Provide one run ID. By default, the driver scans `<base-dir>/UMI/<run-id>/*_S*_L001/input/` for paired FASTQs, determines the library ID from the sample directory prefix, and creates one SSC/DSC job per sample that still needs outputs.

For example, sample directory `18_S1_L001` has library ID `18`.

## Options

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--run-id` | Yes | none | Run identifier used to locate inputs under `<base-dir>/UMI/<run-id>` |
| `--sample-id` | No | discover all samples | Optional exact sample directory to process, such as `18_S1_L001` |
| `--library-id` | No | all discovered libraries | Optional library ID filter, such as `18`; may be repeated or comma-separated |
| `--base-dir` | No | `/data/fadenlab/processed` | Processed data base directory |
| `--r1` | No | inferred from run/sample | Optional full path override for R1 FASTQ.gz input; only valid with `--sample-id` |
| `--r2` | No | inferred from run/sample | Optional full path override for R2 FASTQ.gz input; only valid with `--sample-id` |
| `--output-dir` | No | `<base-dir>/UMI/<run-id>/<sample-id>/ssc_dsc` | Analysis output directory for one explicit sample; only valid with `--sample-id` |
| `--ssc-min-family-size` | No | `5` | SSC family-size threshold for the additional filtered SSC BAM |
| `--ssc-family-size-tag` | No | `cD` | Integer SAM tag used as the SSC family-size value |
| `--reference-genome` | No | `/data/fadenlab/reference/HPV-DeepSeek/HPV_hg38/HPV_hg38.fasta` | HPV_hg38 reference FASTA |
| `--fgbio-jar` | No | `/data/fadenlab/bin/fgbio/fgbio-2.1.0.jar` | fgbio jar path |
| `--umi-extract-script` | No | `extract_umis_by_species.py` next to `ssc_and_dsc.py` | Path to the species-aware UMI extraction script |
| `--umi-catalog-file` | No | `doc/KAPA_Universal_UMI_sequences.txt` | Path to the KAPA Universal UMI catalog file |
| `--three-nt-trim-length` | No | `4` | Total non-genomic prefix length (UMI + T-overhang) to trim for the 3-nt KAPA UMI species |
| `--skip-umi-whitelisting` | No | false | Do not run the index-hopping detection whitelisting step ([details](#index-hopping-detection-umi-whitelisting)) |
| `--umi-whitelist-script` | No | `whitelist_index_hopping_umis.py` next to `ssc_and_dsc.py` | Path to the UMI whitelisting script |
| `--umi-whitelist-max-edit-distance` | No | `1` | Maximum Levenshtein distance accepted as a `near` catalog match |
| `--umi-whitelist-contig-filter` | No | `(?i)hpv` | Regex restricting whitelisting to matching contigs; empty string disables filtering |
| `--skip-index-hopping-audit` | No | false | Do not run the cross-sample index-hopping audit ([details](#index-hopping-detection-cross-sample-audit)); implied by `--skip-umi-whitelisting` |
| `--force-index-hopping-audit` | No | false | Resubmit the audit job even if nothing new was submitted and a prior report exists |
| `--index-hopping-audit-script` | No | `audit_index_hopping.py` next to `ssc_and_dsc.py` | Path to the cross-sample audit script |
| `--index-hopping-min-floor` | No | `10` | Minimum genotype read total below which a non-source sample's reads are always blacklisted |
| `--index-hopping-snr-threshold` | No | `0.01` | Minimum ratio of a non-source sample's genotype total to the source sample's total to rescue it |
| `--index-hopping-cpus` | No | `4` | CPUs for the index-hopping audit SLURM job |
| `--index-hopping-mem` | No | `16G` | SLURM memory for the index-hopping audit job |
| `--index-hopping-time` | No | `04:00:00` | SLURM time limit for the index-hopping audit job |
| `--cpus` | No | `16` | CPUs for BWA and the generated SLURM job |
| `--java-memory` | No | `-Xmx16g` | Java JVM memory option |
| `--partition` | No | `normal` | SLURM partition |
| `--time` | No | `24:00:00` | SLURM time limit |
| `--mem` | No | `64G` | SLURM memory request |
| `--force` | No | false | Create and submit jobs even when expected SSC/DSC outputs already exist |
| `--dry-run` | No | false | Create the job script and log, but do not submit with `sbatch` |
| `--skip-input-check` | No | false | Do not require FASTQ/reference/fgbio paths to exist when generating the job script |

## Default paths

- Base directory: `/data/fadenlab/processed`
- Inferred input directory: `/data/fadenlab/processed/UMI/<run-id>/<sample-id>/input`
- Inferred output directory: `/data/fadenlab/processed/UMI/<run-id>/<sample-id>/ssc_dsc`
- Reference genome: `/data/fadenlab/reference/HPV-DeepSeek/HPV_hg38/HPV_hg38.fasta`
- fgbio jar: `/data/fadenlab/bin/fgbio/fgbio-2.1.0.jar`

## Interactive environment

For interactive testing on ERIS Nucleus:

```bash
salloc --partition=interactive --time=00:30:00 --cpus-per-task=16 --mem=64G
module load Miniforge3/24.11.3-0
source "$EBROOTMINIFORGE3/etc/profile.d/conda.sh"
conda activate fastp-nucleus
```

## UMI species classification

The KAPA Universal UMI Adapter pools 16 UMI sequences of two different lengths (8 that
are 3 nt, 8 that are 5 nt), each followed by a 1-nt T-overhang, allocated independently to
each fragment end. A single fixed `fgbio ExtractUmisFromBam` read structure can only be
exactly correct for one of the two lengths — using one blind guess for every read silently
mis-trims the fragment terminus on whichever species it doesn't match.

Before running `ExtractUmisFromBam`, the pipeline now runs `extract_umis_by_species.py`,
which classifies each mate's UMI species directly from its own leading bases (matched
against `doc/KAPA_Universal_UMI_sequences.txt`, Levenshtein distance ≤ 1 for sequencing-error
tolerance), splits the unmapped BAM into per-species-combination buckets, runs
`ExtractUmisFromBam` once per bucket with the read structure that's exactly correct for that
combination, and merges the results back into a single UMI-tagged unmapped BAM before the
rest of the pipeline continues unchanged. Read pairs where either mate's UMI doesn't match
any known sequence within tolerance are routed to a separate `unclassified` bucket and
excluded, rather than guessed.

`--three-nt-trim-length` (default `4`) controls the total non-genomic prefix trimmed for the
3-nt species (UMI + 1-base T-overhang, no padding assumed); it's exposed as a parameter
because this default has not been empirically validated against real sequencing data.
Requires `pysam` and `python-Levenshtein` in the `fastp-nucleus` conda environment.

## Index-hopping detection: UMI whitelisting

After the SSC and DSC BAMs are finalized, the pipeline runs `whitelist_index_hopping_umis.py`
on each one: for every template mapped to an HPV contig, it classifies the `RX:Z:` UMI tag
against the KAPA Universal UMI catalog using the same Levenshtein-distance matcher
`extract_umis_by_species.py` uses before extraction (shared in `umi_matching.py`), rather than
a Hamming/positional-mismatch whitelist. Levenshtein distance is length-aware (each observed
UMI segment is checked only against the catalog subset of the same length) and correctly
scores insertion/deletion errors that a positional mismatch count cannot detect at all.

Each read is classified as `exact` (UMI matches the catalog exactly), `near` (within
`--umi-whitelist-max-edit-distance`, default 1), `garbage` (further than that), or
`unclassified_length` (the observed UMI segment is neither 3 nor 5 nt). The result is a
per-read report — `<sample_id>.SSC_umi_whitelist.tsv` / `<sample_id>.DSC_umi_whitelist.tsv` in
`output_SSC/` / `output_DSC/` — with columns for the genomic signature (contig, position, mate
position, template length) and the UMI match. This is the whitelisting step only: it reports
per-read UMI/signature classification and does not itself compare across samples, apply any
cross-sample accept/reject threshold, or modify the input BAM.

Options:

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--skip-umi-whitelisting` | No | false | Do not run the whitelisting step |
| `--umi-whitelist-script` | No | `whitelist_index_hopping_umis.py` next to `ssc_and_dsc.py` | Path to the whitelisting script |
| `--umi-whitelist-max-edit-distance` | No | `1` | Maximum Levenshtein distance accepted as a `near` catalog match |
| `--umi-whitelist-contig-filter` | No | `(?i)hpv` | Regex restricting whitelisting to matching contigs; empty string disables filtering |

## Index-hopping detection: cross-sample audit

After all of a run's SSC/DSC sample jobs are queued, the pipeline submits one additional
SLURM job -- with `--dependency=afterany:<job IDs>` on the samples this invocation just
submitted -- that runs `audit_index_hopping.py` once for SSC and once for DSC. It reads every
sample's UMI whitelist report (see above) for the run and groups `exact`/`near` reads across
*all* samples by physical signature (contig, position, mate position, template length) plus
canonical UMI. Two different samples producing a read with an identical signature and UMI is,
for practical purposes, the same physical molecule read under two different sample indices --
index hopping -- rather than coincidence.

For each such cross-sample group, the sample with the most total UMI-matched reads for that
contig is treated as the true source. Every other sample's reads in the group are then
audited: below `--index-hopping-min-floor` total reads for that contig, or below
`--index-hopping-snr-threshold` of the source sample's total, they are blacklisted as
index-hopping artifacts; otherwise they are rescued as a plausible low-level co-infection.
Blacklisted reads (both mates) are stripped to produce `*_clean.bam` files -- for SSC, both
the unfiltered and family-size-filtered BAMs get a clean variant, using the same blacklist --
alongside a per-run `IHOP_Final_Report.SSC.tsv` / `IHOP_Final_Report.DSC.tsv` under
`<run-id>/index_hopping/results/`.

Known constraints, by design:

- **Whole-run manifest, per-invocation dependency.** The manifest always covers every sample
  discovered in the run directory, regardless of any `--library-id` filter used for this
  invocation -- but the SLURM dependency only waits on jobs *this invocation* submitted.
  Samples from a different library submitted in an earlier or later invocation are included
  in the manifest by their expected path; if that path doesn't exist yet when the audit job
  actually runs, that sample is skipped with a warning rather than failing the audit. Re-run
  once every library has been submitted to get a complete audit.
- **Whole-run only.** A `--sample-id`-scoped invocation never triggers or updates the audit
  (there's no "rest of the run" to compare against, and `--output-dir` overrides mean its
  outputs may not even be at the path a whole-run scan would find). Re-run without
  `--sample-id` to refresh the audit after a targeted single-sample re-run.
- **Idempotent by default.** If nothing new was submitted this invocation and a prior report
  already exists, the audit job is not resubmitted; pass `--force-index-hopping-audit` (or
  `--force`) to force it anyway.
- **DSC audited independently of SSC.** `GroupReadsByUmi`'s `adjacency` (SSC) and `paired`
  (DSC) strategies produce unrelated qname spaces, so SSC and DSC signatures/read names are
  never joined against each other.

Options:

| Option | Required | Default | Description |
|--------|----------|---------|-------------|
| `--skip-index-hopping-audit` | No | false | Do not run the cross-sample audit (implied by `--skip-umi-whitelisting`) |
| `--force-index-hopping-audit` | No | false | Resubmit the audit job even if nothing new was submitted and a prior report exists |
| `--index-hopping-audit-script` | No | `audit_index_hopping.py` next to `ssc_and_dsc.py` | Path to the audit script |
| `--index-hopping-min-floor` | No | `10` | Minimum genotype read total below which a non-source sample's reads are always blacklisted |
| `--index-hopping-snr-threshold` | No | `0.01` | Minimum ratio of a non-source sample's genotype total to the source sample's total to rescue it |
| `--index-hopping-cpus` / `--index-hopping-mem` / `--index-hopping-time` | No | `4` / `16G` / `04:00:00` | SLURM resources for the audit job |

## SSC family-size filtering

The Python port keeps the original unfiltered SSC BAM and adds a second SSC BAM filtered by consensus family size.

Default behavior:

```bash
--ssc-min-family-size 5 --ssc-family-size-tag cD
```

The generated SLURM job keeps reads only when the selected integer SAM tag is greater than the threshold. With the defaults, the extra filtered SSC BAM keeps reads with `cD > 5`.

The underscore aliases `--run_id`, `--library_id`, `--base_dir`, and `--output_dir` are also accepted.

## Example run dry run

```bash
cd /data/fadenlab/scripts/scripts-nucleus/reads_coverage

python ssc_and_dsc.py \
    --run-id 10_27_2025-test \
    --ssc-min-family-size 5 \
    --dry-run
```

For the example folder shown below, the script detects paired FASTQs for `18_S1_L001`, `19_S1_L001`, and `4_S1_L001`. Because `4_S1_L001` already has expected SSC/DSC outputs, it is skipped unless `--force` is used, so the dry run creates jobs for `18_S1_L001` and `19_S1_L001`. Because this is a whole-run invocation (no `--sample-id`) with 3 samples discovered, it also creates an [index-hopping audit](#index-hopping-detection-cross-sample-audit) manifest and job script covering all 3 samples -- including `4_S1_L001`, even though no new job was created for it this run.

```text
/data/fadenlab/processed/UMI/10_27_2025-test/18_S1_L001/input/18_S1_L001_R1_001.fastq.gz
/data/fadenlab/processed/UMI/10_27_2025-test/18_S1_L001/input/18_S1_L001_R2_001.fastq.gz
```

Dry run creates the job scripts and master logs, but does not call `sbatch`.

## Filtering by library ID

To test one library only:

```bash
python ssc_and_dsc.py \
    --run-id 10_27_2025-test \
    --library-id 18 \
    --dry-run
```

`--library-id` only filters which samples get an SSC/DSC job *this invocation* -- the
[index-hopping audit](#index-hopping-detection-cross-sample-audit)'s manifest always covers
every sample already discoverable in the run directory, regardless of this filter.

To process one exact sample directory:

```bash
python ssc_and_dsc.py \
    --run-id 10_27_2025-test \
    --sample-id 18_S1_L001 \
    --dry-run
```

## Output layout

The output directory contains:

- `job_scripts/`: generated SLURM job script
- `logs/`: master analysis log plus SLURM stdout/stderr
- `qc/`: fastp JSON, HTML, and log files
- `output_SSC/`: unfiltered SSC BAM/index/coverage, family-size-filtered SSC BAM/index/coverage, and the SSC UMI whitelist report
- `output_DSC/`: DSC BAM/index/coverage, duplex metrics, and the DSC UMI whitelist report

At the run level (`<base-dir>/UMI/<run-id>/index_hopping/`), the cross-sample audit adds:

- `sample_bam_manifest.tsv`: whole-run sample/BAM/whitelist-report manifest
- `job_scripts/`, `logs/`: the audit SLURM job script and its stdout/stderr
- `results/`: `*_clean.bam` (+ index) per sample/BAM-tier and `IHOP_Final_Report.SSC.tsv` / `IHOP_Final_Report.DSC.tsv`

The master log records full paths for all inputs and expected outputs, all run parameters, the ERIS Nucleus module names, and expected program versions. Runtime tool version output is appended by the SLURM job.

## Nucleus software stack

The per-sample SSC/DSC job loads these Nucleus modules:

- `Java/17.0.15`
- `GATK/4.6.1.0-GCCcore-13.3.0-Java-17`
- `BWA/0.7.19-GCCcore-13.3.0`
- `SAMtools/1.22.1-GCC-13.3.0`
- `BEDTools/2.31.1-GCC-13.3.0`
- `Miniforge3/24.11.3-0`

It then activates the `fastp-nucleus` conda environment, which also provides `pysam` and
`python-Levenshtein` for the UMI species classification and whitelisting steps.

The index-hopping cross-sample audit job (`audit_index_hopping.py`) only needs `pysam` (BAM
filtering uses pysam's bundled htslib, not the system `samtools`), so it loads only
`Miniforge3/24.11.3-0` and activates `fastp-nucleus` -- it does not load GATK, BWA, SAMtools,
or BEDTools.
