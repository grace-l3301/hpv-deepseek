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

For the example folder shown below, the script detects paired FASTQs for `18_S1_L001`, `19_S1_L001`, and `4_S1_L001`. Because `4_S1_L001` already has expected SSC/DSC outputs, it is skipped unless `--force` is used, so the dry run creates jobs for `18_S1_L001` and `19_S1_L001`.

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
- `output_SSC/`: unfiltered SSC BAM/index/coverage and family-size-filtered SSC BAM/index/coverage
- `output_DSC/`: DSC BAM/index/coverage and duplex metrics

The master log records full paths for all inputs and expected outputs, all run parameters, the ERIS Nucleus module names, and expected program versions. Runtime tool version output is appended by the SLURM job.

## Nucleus software stack

The generated job loads these Nucleus modules:

- `Java/17.0.15`
- `GATK/4.6.1.0-GCCcore-13.3.0-Java-17`
- `BWA/0.7.19-GCCcore-13.3.0`
- `SAMtools/1.22.1-GCC-13.3.0`
- `BEDTools/2.31.1-GCC-13.3.0`
- `Miniforge3/24.11.3-0`

It then activates the `fastp-nucleus` conda environment.
