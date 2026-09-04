# reads_coverage.py

Generate and optionally submit ERIS Nucleus SLURM jobs for the UMI deduplication and coverage pipeline.

**Version:** 9.1.0

## Overview

`reads_coverage.py` scans one processed UMI run directory, detects sample directories, writes one SLURM job script per sample, and submits those jobs unless `--dry-run` is used.

This is the batch UMI pipeline that produces the per-sample BAM files used by downstream scripts:

```text
/data/fadenlab/processed/UMI/<RUN_DATE>/<sample>/output/<sample>_umi_dedup_sorted.bam
```

## Usage

The run date is currently configured inside the script:

```python
RUN_DATE = "10_27_2025-test"
```

Update `RUN_DATE` in `reads_coverage.py`, then run from the Nucleus scripts directory:

```bash
cd /data/fadenlab/scripts/scripts-nucleus/reads_coverage

# Create job scripts and logs without submitting jobs
python reads_coverage.py --dry-run

# Create and submit one SLURM job per detected sample
python reads_coverage.py
```

## Arguments

| Option | Description |
|--------|-------------|
| `--dry-run` | Create job scripts but skip `sbatch` submission |

## Hardcoded defaults

| Variable | Value |
|----------|-------|
| `BASE_DIR` | `/data/fadenlab/processed` |
| `RUN_DATE` | `10_27_2025-test` |
| `REFERENCE_GENOME` | `/data/fadenlab/reference/HPV-DeepSeek/HPV_hg38/HPV_hg38.fasta` |
| `FGBIO_JAR` | `/data/fadenlab/bin/fgbio/fgbio-2.1.0.jar` |
| `THREAD_COUNT` | `2` |
| `JVM_OPTIONS` | `-Xmx32g` |

## Input requirements

The script scans:

```text
/data/fadenlab/processed/UMI/<RUN_DATE>/*_S*_L001
```

Detected sample directories must match a numeric sample pattern such as `10_S41_L001` and contain paired FASTQ inputs:

```text
/data/fadenlab/processed/UMI/<RUN_DATE>/<sample>/input/<sample>_R1_001.fastq.gz
/data/fadenlab/processed/UMI/<RUN_DATE>/<sample>/input/<sample>_R2_001.fastq.gz
```

Each sample directory should also contain an `output/` directory before submission because the generated SLURM script writes stdout, stderr, fastp reports, BAMs, and coverage files there.

## Generated files

The driver creates run-level files under:

```text
/data/fadenlab/processed/UMI/<RUN_DATE>/
```

| Path | Description |
|------|-------------|
| `job_scripts/job_<sample>.sh` | Generated per-sample SLURM script |
| `pipeline_log.txt` | Run-level log with detected samples and expected versions |

Each per-sample job writes outputs under:

```text
/data/fadenlab/processed/UMI/<RUN_DATE>/<sample>/output/
```

| Output | Description |
|--------|-------------|
| `<sample>_umi_dedup_sorted.bam` | UMI consensus deduplicated BAM |
| `<sample>_umi_dedup_sorted.bam.bai` | BAM index |
| `<sample>.coverage.txt` | Coverage summary from the unfiltered BAM |
| `<sample>.filtered.coverage.txt` | Coverage summary after removing `HPV31_Ref:4093-4137` |
| `<sample>_fastp.html` | fastp HTML report |
| `<sample>_fastp.log` | fastp log |
| `umi_pipeline_<jobid>.out` | SLURM stdout |
| `umi_pipeline_<jobid>.err` | SLURM stderr |

## SLURM configuration

Each generated per-sample job uses:

| Resource | Value |
|----------|-------|
| CPUs | 13 |
| Memory | 62G |
| Time | 24 hours |
| Partition | cluster default unless overridden externally |

## Nucleus software stack

Generated jobs load these ERIS Nucleus modules:

- `GATK/4.6.1.0-GCCcore-13.3.0-Java-17`
- `BWA/0.7.19-GCCcore-13.3.0`
- `SAMtools/1.22.1-GCC-13.3.0`
- `BEDTools/2.31.1-GCC-13.3.0`
- `Miniforge3/24.11.3-0`

The job template initializes conda with:

```bash
source "$EBROOTMINIFORGE3/etc/profile.d/conda.sh"
```

## Conda environment

Generated jobs activate the Nucleus fastp environment:

```bash
conda activate fastp-nucleus
```

The Nucleus environment file is `envs/fastp.yml`, and the environment name defined inside that file is `fastp-nucleus`.

## Workflow

For each detected sample, the generated SLURM job:

1. Converts paired FASTQ inputs to an unmapped BAM with GATK `FastqToSam`
2. Extracts UMIs with fgbio `ExtractUmisFromBam`
3. Converts UMI-extracted BAM back to FASTQ
4. Runs `fastp` adapter trimming and quality filtering
5. Aligns reads to HPV_hg38 with BWA
6. Merges aligned reads with UMI metadata
7. Groups reads by UMI with fgbio `GroupReadsByUmi`
8. Calls molecular consensus reads with fgbio `CallMolecularConsensusReads`
9. Realigns consensus reads and writes `<sample>_umi_dedup_sorted.bam`
10. Generates coverage summaries before and after the HPV31 exclusion filter