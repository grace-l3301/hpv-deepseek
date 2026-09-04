#!/usr/bin/env python3
"""
ssc_and_dsc.py

Generate and optionally submit Nucleus-compatible SLURM jobs for SSC and DSC
UMI consensus analysis from paired-end FASTQ samples in one run folder.

Author: Samuli Eldfors
"""

__version__ = "1.4.1"

import argparse
import re
import shlex
import stat
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from string import Template


DEFAULT_REFERENCE_GENOME = "/data/fadenlab/reference/HPV-DeepSeek/HPV_hg38/HPV_hg38.fasta"
DEFAULT_FGBIO_JAR = "/data/fadenlab/bin/fgbio/fgbio-2.1.0.jar"
DEFAULT_BASE_DIR = "/data/fadenlab/processed"
DEFAULT_CPUS = 16
DEFAULT_JVM_OPTIONS = "-Xmx16g"
DEFAULT_SLURM_MEM = "64G"
DEFAULT_SLURM_TIME = "24:00:00"
DEFAULT_SLURM_PARTITION = "normal"
DEFAULT_SSC_MIN_FAMILY_SIZE = 5
DEFAULT_SSC_FAMILY_SIZE_TAG = "cD"
DEFAULT_UMI_EXTRACT_SCRIPT = str(Path(__file__).resolve().parent / "extract_umis_by_species.py")
DEFAULT_UMI_CATALOG_FILE = str(Path(__file__).resolve().parent / "doc" / "KAPA_Universal_UMI_sequences.txt")
DEFAULT_THREE_NT_TRIM_LENGTH = 4

NUCLEUS_MODULES = {
    "Java": "Java/17.0.15",
    "GATK": "GATK/4.6.1.0-GCCcore-13.3.0-Java-17",
    "BWA": "BWA/0.7.19-GCCcore-13.3.0",
    "SAMtools": "SAMtools/1.22.1-GCC-13.3.0",
    "BEDTools": "BEDTools/2.31.1-GCC-13.3.0",
    "Miniforge3": "Miniforge3/24.11.3-0",
}

PROGRAM_VERSIONS = {
    "Java": "17.0.15",
    "GATK": "4.6.1.0",
    "BWA": "0.7.19",
    "SAMtools": "1.22.1",
    "BEDTools": "2.31.1",
    "fgbio": "2.1.0",
    "fastp": "from conda environment fastp-nucleus",
    "pysam": "from conda environment fastp-nucleus",
    "python-Levenshtein": "from conda environment fastp-nucleus",
}

SAMPLE_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
LIBRARY_ID_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")
SAMPLE_DIR_PATTERN = re.compile(r"^(?P<library_id>[^_]+)_S\d+_L001$")
R1_FASTQ_PATTERN = re.compile(r"^(?P<sample_id>.+)_R1_001\.fastq\.gz$")


@dataclass(frozen=True)
class SampleJob:
    sample_id: str
    library_id: str
    r1_fastq: Path
    r2_fastq: Path


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate and optionally submit Nucleus SLURM jobs for SSC/DSC "
            "UMI consensus analysis samples in one run folder."
        )
    )
    parser.add_argument(
        "--run-id",
        "--run_id",
        dest="run_id",
        required=True,
        help="Run identifier used to locate inputs under <base-dir>/UMI/<run-id>.",
    )
    parser.add_argument(
        "--sample-id",
        help=(
            "Optional exact sample directory to process, for example 18_S1_L001. "
            "If omitted, all valid samples in the run folder are discovered."
        ),
    )
    parser.add_argument(
        "--library-id",
        "--library_id",
        dest="library_ids",
        action="append",
        default=[],
        help=(
            "Optional library ID filter, for example 18. May be repeated or given as "
            "a comma-separated list. If omitted, all discovered library IDs are included."
        ),
    )
    parser.add_argument(
        "--base-dir",
        "--base_dir",
        dest="base_dir",
        default=DEFAULT_BASE_DIR,
        help=f"Processed data base directory. Default: {DEFAULT_BASE_DIR}.",
    )
    parser.add_argument(
        "--r1",
        help=(
            "Optional full path override for R1 FASTQ.gz. By default this is inferred as "
            "<base-dir>/UMI/<run-id>/<sample-id>/input/<sample-id>_R1_001.fastq.gz. "
            "Only valid with --sample-id."
        ),
    )
    parser.add_argument(
        "--r2",
        help=(
            "Optional full path override for R2 FASTQ.gz. By default this is inferred as "
            "<base-dir>/UMI/<run-id>/<sample-id>/input/<sample-id>_R2_001.fastq.gz. "
            "Only valid with --sample-id."
        ),
    )
    parser.add_argument(
        "--output-dir",
        "--output_dir",
        dest="output_dir",
        help=(
            "Analysis output directory for this sample. Subdirectories are created inside it. "
            "Default: <base-dir>/UMI/<run-id>/<sample-id>/ssc_dsc. Only valid with --sample-id."
        ),
    )
    parser.add_argument(
        "--ssc-min-family-size",
        type=int,
        default=DEFAULT_SSC_MIN_FAMILY_SIZE,
        help=(
            "SSC family-size threshold for the additional filtered SSC BAM. "
            "Reads are retained only when the selected family-size tag is greater than this value. "
            f"Default: {DEFAULT_SSC_MIN_FAMILY_SIZE}."
        ),
    )
    parser.add_argument(
        "--ssc-family-size-tag",
        default=DEFAULT_SSC_FAMILY_SIZE_TAG,
        help=(
            "Integer SAM tag used as the SSC family-size value in the final SSC BAM. "
            f"Default: {DEFAULT_SSC_FAMILY_SIZE_TAG}."
        ),
    )
    parser.add_argument("--reference-genome", default=DEFAULT_REFERENCE_GENOME, help="Full path to HPV_hg38 FASTA.")
    parser.add_argument("--fgbio-jar", default=DEFAULT_FGBIO_JAR, help="Full path to the fgbio jar.")
    parser.add_argument(
        "--umi-extract-script",
        default=DEFAULT_UMI_EXTRACT_SCRIPT,
        help="Full path to the species-aware UMI extraction script (extract_umis_by_species.py).",
    )
    parser.add_argument(
        "--umi-catalog-file",
        default=DEFAULT_UMI_CATALOG_FILE,
        help="Full path to the KAPA Universal UMI catalog file (bare UMI + full-sequence columns).",
    )
    parser.add_argument(
        "--three-nt-trim-length",
        type=int,
        default=DEFAULT_THREE_NT_TRIM_LENGTH,
        help=(
            "Total non-genomic prefix length (UMI + T-overhang) to trim for the 3-nt KAPA UMI "
            f"species. Default: {DEFAULT_THREE_NT_TRIM_LENGTH} (3-nt UMI + 1-base T-overhang, no padding)."
        ),
    )
    parser.add_argument("--cpus", type=int, default=DEFAULT_CPUS, help=f"CPUs for BWA and SLURM. Default: {DEFAULT_CPUS}.")
    parser.add_argument("--java-memory", default=DEFAULT_JVM_OPTIONS, help=f"Java JVM memory option. Default: {DEFAULT_JVM_OPTIONS}.")
    parser.add_argument("--partition", default=DEFAULT_SLURM_PARTITION, help=f"SLURM partition. Default: {DEFAULT_SLURM_PARTITION}.")
    parser.add_argument("--time", default=DEFAULT_SLURM_TIME, help=f"SLURM time limit. Default: {DEFAULT_SLURM_TIME}.")
    parser.add_argument("--mem", default=DEFAULT_SLURM_MEM, help=f"SLURM memory. Default: {DEFAULT_SLURM_MEM}.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Create and submit jobs even when expected SSC/DSC output files already exist.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Create job script and log, but do not submit with sbatch.")
    parser.add_argument(
        "--skip-input-check",
        action="store_true",
        help="Do not require FASTQ/reference/fgbio paths to exist when generating the job script.",
    )
    return parser.parse_args()


def resolve_path(path_text):
    return str(Path(path_text).expanduser().resolve())


def normalize_library_ids(library_ids):
    normalized = []
    for library_id_text in library_ids or []:
        normalized.extend(library_id.strip() for library_id in library_id_text.split(",") if library_id.strip())
    return normalized


def extract_library_id(sample_id):
    match = SAMPLE_DIR_PATTERN.match(sample_id)
    if not match:
        return None
    return match.group("library_id")


def build_run_dir(args):
    base_dir = Path(resolve_path(args.base_dir))
    return base_dir / "UMI" / args.run_id


def build_run_paths(args, sample_id=None, r1_fastq=None, r2_fastq=None):
    sample_id = sample_id or args.sample_id
    if not sample_id:
        raise ValueError("sample_id is required to build sample paths.")
    base_dir = Path(resolve_path(args.base_dir))
    run_dir = build_run_dir(args)
    sample_dir = run_dir / sample_id
    input_dir = sample_dir / "input"
    r1_fastq = Path(resolve_path(r1_fastq or args.r1)) if r1_fastq or args.r1 else input_dir / f"{sample_id}_R1_001.fastq.gz"
    r2_fastq = Path(resolve_path(r2_fastq or args.r2)) if r2_fastq or args.r2 else input_dir / f"{sample_id}_R2_001.fastq.gz"
    return {
        "base_dir": base_dir,
        "run_dir": run_dir,
        "sample_dir": sample_dir,
        "input_dir": input_dir,
        "r1_fastq": r1_fastq,
        "r2_fastq": r2_fastq,
    }


def validate_args(args):
    if not RUN_ID_PATTERN.match(args.run_id):
        raise ValueError("--run-id may contain only letters, numbers, '.', '_', and '-'.")
    if args.sample_id and not SAMPLE_ID_PATTERN.match(args.sample_id):
        raise ValueError("--sample-id may contain only letters, numbers, '.', '_', and '-'.")
    for library_id in normalize_library_ids(args.library_ids):
        if not LIBRARY_ID_PATTERN.match(library_id):
            raise ValueError("--library-id may contain only letters, numbers, '.', '_', and '-'.")
    if (args.r1 or args.r2) and not args.sample_id:
        raise ValueError("--r1 and --r2 overrides may only be used with --sample-id.")
    if (args.r1 and not args.r2) or (args.r2 and not args.r1):
        raise ValueError("--r1 and --r2 must be provided together.")
    if args.output_dir and not args.sample_id:
        raise ValueError("--output-dir may only be used with --sample-id.")
    if args.ssc_min_family_size < 0:
        raise ValueError("--ssc-min-family-size must be >= 0.")
    if not re.match(r"^[A-Za-z][A-Za-z0-9]$", args.ssc_family_size_tag):
        raise ValueError("--ssc-family-size-tag must look like a two-character SAM tag, for example cD or cM.")
    if args.cpus < 1:
        raise ValueError("--cpus must be >= 1.")
    if args.three_nt_trim_length < 4:
        raise ValueError(
            "--three-nt-trim-length must be >= 4 (3-nt UMI plus at least a 1-base T-overhang skip)."
        )

    run_dir = build_run_dir(args)
    if not args.sample_id or not args.skip_input_check:
        if not run_dir.is_dir():
            raise FileNotFoundError(f"run directory not found: {run_dir}")

    if not args.skip_input_check:
        for label, path_text in (
            ("reference genome", args.reference_genome),
            ("fgbio jar", args.fgbio_jar),
            ("UMI extraction script", args.umi_extract_script),
            ("UMI catalog file", args.umi_catalog_file),
        ):
            if not Path(path_text).expanduser().exists():
                raise FileNotFoundError(f"{label} not found: {path_text}")


def discover_sample_jobs(args):
    selected_library_ids = set(normalize_library_ids(args.library_ids))

    if args.sample_id:
        library_id = extract_library_id(args.sample_id) or args.sample_id
        if selected_library_ids and library_id not in selected_library_ids:
            return [], [f"{args.sample_id}: sample library ID {library_id} does not match --library-id filter"]
        run_paths = build_run_paths(args)
        return [
            SampleJob(
                sample_id=args.sample_id,
                library_id=library_id,
                r1_fastq=run_paths["r1_fastq"],
                r2_fastq=run_paths["r2_fastq"],
            )
        ], []

    run_dir = build_run_dir(args)
    sample_jobs = []
    skipped_inputs = []
    seen_sample_ids = set()

    for r1_fastq in sorted(run_dir.glob("*_S*_L001/input/*_R1_001.fastq.gz")):
        match = R1_FASTQ_PATTERN.match(r1_fastq.name)
        if not match:
            continue

        sample_id = match.group("sample_id")
        if sample_id in seen_sample_ids:
            continue
        seen_sample_ids.add(sample_id)

        sample_dir = r1_fastq.parent.parent
        if sample_dir.name != sample_id:
            skipped_inputs.append(f"{r1_fastq}: sample directory name does not match FASTQ prefix")
            continue

        sample_match = SAMPLE_DIR_PATTERN.match(sample_id)
        if not sample_match:
            skipped_inputs.append(f"{sample_dir}: sample directory does not match <library_id>_S*_L001")
            continue

        library_id = sample_match.group("library_id")
        if selected_library_ids and library_id not in selected_library_ids:
            continue

        r2_fastq = r1_fastq.parent / f"{sample_id}_R2_001.fastq.gz"
        if not r2_fastq.exists():
            skipped_inputs.append(f"{sample_id}: missing R2 FASTQ {r2_fastq}")
            continue

        sample_jobs.append(
            SampleJob(
                sample_id=sample_id,
                library_id=library_id,
                r1_fastq=r1_fastq,
                r2_fastq=r2_fastq,
            )
        )

    return sample_jobs, skipped_inputs


def build_sample_args(args, sample_job):
    sample_args = argparse.Namespace(**vars(args))
    sample_args.sample_id = sample_job.sample_id
    sample_args.r1 = str(sample_job.r1_fastq)
    sample_args.r2 = str(sample_job.r2_fastq)
    if not args.sample_id:
        sample_args.output_dir = None
    return sample_args


def validate_sample_inputs(args, sample_job):
    if args.skip_input_check:
        return
    input_dir = sample_job.r1_fastq.parent
    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {input_dir}")
    for label, fastq_path in (("R1 FASTQ", sample_job.r1_fastq), ("R2 FASTQ", sample_job.r2_fastq)):
        if not fastq_path.exists():
            raise FileNotFoundError(f"{label} not found: {fastq_path}")


def build_paths(args):
    sample_id = args.sample_id
    run_paths = build_run_paths(args)
    analysis_dir = Path(resolve_path(args.output_dir)) if args.output_dir else run_paths["sample_dir"] / "ssc_dsc"
    paths = {
        "base_dir": run_paths["base_dir"],
        "run_dir": run_paths["run_dir"],
        "sample_dir": run_paths["sample_dir"],
        "input_dir": run_paths["input_dir"],
        "r1_fastq": run_paths["r1_fastq"],
        "r2_fastq": run_paths["r2_fastq"],
        "analysis_dir": analysis_dir,
        "job_dir": analysis_dir / "job_scripts",
        "log_dir": analysis_dir / "logs",
        "qc_dir": analysis_dir / "qc",
        "temp_dir": analysis_dir / "temp",
        "temp_ssc_dir": analysis_dir / "temp_SSC",
        "temp_dsc_dir": analysis_dir / "temp_DSC",
        "output_ssc_dir": analysis_dir / "output_SSC",
        "output_dsc_dir": analysis_dir / "output_DSC",
    }
    paths["job_script"] = paths["job_dir"] / f"job_{sample_id}_ssc_dsc.sh"
    paths["analysis_log"] = paths["log_dir"] / f"{sample_id}_ssc_dsc_analysis_log.txt"
    paths["slurm_stdout"] = paths["log_dir"] / f"{sample_id}_ssc_dsc_%j.out"
    paths["slurm_stderr"] = paths["log_dir"] / f"{sample_id}_ssc_dsc_%j.err"
    paths["fastp_json"] = paths["qc_dir"] / f"{sample_id}_fastp.json"
    paths["fastp_html"] = paths["qc_dir"] / f"{sample_id}_fastp.html"
    paths["fastp_log"] = paths["qc_dir"] / f"{sample_id}_fastp.log"
    paths["ssc_bam"] = paths["output_ssc_dir"] / f"{sample_id}_umi_dedup_sorted_SSC.bam"
    paths["ssc_bai"] = Path(str(paths["ssc_bam"]) + ".bai")
    paths["ssc_coverage"] = paths["output_ssc_dir"] / f"{sample_id}.SSC_coverage.txt"
    paths["ssc_filtered_bam"] = paths["output_ssc_dir"] / (
        f"{sample_id}_umi_dedup_sorted_SSC.{args.ssc_family_size_tag}_gt_{args.ssc_min_family_size}.bam"
    )
    paths["ssc_filtered_bai"] = Path(str(paths["ssc_filtered_bam"]) + ".bai")
    paths["ssc_filtered_coverage"] = paths["output_ssc_dir"] / (
        f"{sample_id}.SSC_{args.ssc_family_size_tag}_gt_{args.ssc_min_family_size}_coverage.txt"
    )
    paths["dsc_metrics_prefix"] = paths["output_dsc_dir"] / sample_id
    paths["dsc_bam"] = paths["output_dsc_dir"] / f"{sample_id}_umi_dedup_sorted_DSC.bam"
    paths["dsc_bai"] = Path(str(paths["dsc_bam"]) + ".bai")
    paths["dsc_coverage"] = paths["output_dsc_dir"] / f"{sample_id}.DSC_coverage.txt"
    return paths


def expected_completion_files(paths):
    return (
        paths["ssc_bam"],
        paths["ssc_bai"],
        paths["ssc_coverage"],
        paths["ssc_filtered_bam"],
        paths["ssc_filtered_bai"],
        paths["ssc_filtered_coverage"],
        paths["dsc_bam"],
        paths["dsc_bai"],
        paths["dsc_coverage"],
    )


def outputs_are_complete(paths):
    return all(path.is_file() and path.stat().st_size > 0 for path in expected_completion_files(paths))


def create_directories(paths):
    for key in ("analysis_dir", "job_dir", "log_dir", "qc_dir", "output_ssc_dir", "output_dsc_dir"):
        paths[key].mkdir(parents=True, exist_ok=True)


def write_analysis_log(args, paths):
    reference_genome = resolve_path(args.reference_genome)
    fgbio_jar = resolve_path(args.fgbio_jar)
    umi_extract_script = resolve_path(args.umi_extract_script)
    umi_catalog_file = resolve_path(args.umi_catalog_file)

    with paths["analysis_log"].open("w") as log_file:
        log_file.write(f"Script name: ssc_and_dsc-v{__version__}.py\n")
        log_file.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write("Cluster: ERIS Nucleus\n")
        log_file.write(f"Base directory: {paths['base_dir']}\n")
        log_file.write(f"Run ID: {args.run_id}\n")
        log_file.write(f"Run directory: {paths['run_dir']}\n")
        log_file.write(f"Sample ID: {args.sample_id}\n")
        log_file.write(f"Sample directory: {paths['sample_dir']}\n")
        log_file.write(f"Input directory: {paths['input_dir']}\n")
        log_file.write(f"Analysis directory: {paths['analysis_dir']}\n")
        log_file.write(f"Dry run: {args.dry_run}\n")
        log_file.write("\nParameters:\n")
        log_file.write(f"  CPUs: {args.cpus}\n")
        log_file.write(f"  Java memory: {args.java_memory}\n")
        log_file.write(f"  SLURM partition: {args.partition}\n")
        log_file.write(f"  SLURM time: {args.time}\n")
        log_file.write(f"  SLURM memory: {args.mem}\n")
        log_file.write(f"  SSC family-size tag: {args.ssc_family_size_tag}\n")
        log_file.write(f"  SSC family-size filter: {args.ssc_family_size_tag} > {args.ssc_min_family_size}\n")
        log_file.write(f"  3-nt UMI species trim length: {args.three_nt_trim_length}\n")
        log_file.write("\nInput files:\n")
        log_file.write(f"  R1 FASTQ: {paths['r1_fastq']}\n")
        log_file.write(f"  R2 FASTQ: {paths['r2_fastq']}\n")
        log_file.write(f"  Reference genome: {reference_genome}\n")
        log_file.write(f"  fgbio jar: {fgbio_jar}\n")
        log_file.write(f"  UMI extraction script: {umi_extract_script}\n")
        log_file.write(f"  UMI catalog file: {umi_catalog_file}\n")
        log_file.write("\nGenerated files:\n")
        for label in (
            "job_script",
            "analysis_log",
            "slurm_stdout",
            "slurm_stderr",
            "fastp_json",
            "fastp_html",
            "fastp_log",
            "ssc_bam",
            "ssc_bai",
            "ssc_coverage",
            "ssc_filtered_bam",
            "ssc_filtered_bai",
            "ssc_filtered_coverage",
            "dsc_metrics_prefix",
            "dsc_bam",
            "dsc_bai",
            "dsc_coverage",
        ):
            log_file.write(f"  {label}: {paths[label]}\n")
        log_file.write("\nNucleus modules:\n")
        for tool, module_name in NUCLEUS_MODULES.items():
            log_file.write(f"  {tool}: {module_name}\n")
        log_file.write("\nProgram versions:\n")
        for tool, version in PROGRAM_VERSIONS.items():
            log_file.write(f"  {tool}: {version}\n")
        log_file.write("\nRuntime version output will be appended by the SLURM job.\n")


def shell_quote(path_or_value):
    return shlex.quote(str(path_or_value))


def build_job_script(args, paths):
    template = Template(r'''#!/bin/bash
#SBATCH --job-name=ssc_dsc_${sample_id}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=${cpus}
#SBATCH --mem=${slurm_mem}
#SBATCH --time=${slurm_time}
#SBATCH --partition=${partition}
#SBATCH --output=${slurm_stdout}
#SBATCH --error=${slurm_stderr}

set -euo pipefail

SAMPLE_ID=${sample_id_quoted}
R1_FASTQ=${r1_fastq}
R2_FASTQ=${r2_fastq}
REFERENCE_GENOME=${reference_genome}
FGBIO_JAR=${fgbio_jar}
UMI_EXTRACT_SCRIPT=${umi_extract_script}
UMI_CATALOG_FILE=${umi_catalog_file}
THREE_NT_TRIM_LENGTH=${three_nt_trim_length}
JVM_OPTIONS=${jvm_options}
CPU=${cpus}
ANALYSIS_DIR=${analysis_dir}
TEMP_DIR=${temp_dir}
TEMP_SSC_DIR=${temp_ssc_dir}
TEMP_DSC_DIR=${temp_dsc_dir}
OUTPUT_SSC_DIR=${output_ssc_dir}
OUTPUT_DSC_DIR=${output_dsc_dir}
QC_DIR=${qc_dir}
ANALYSIS_LOG=${analysis_log}
FASTP_JSON=${fastp_json}
FASTP_HTML=${fastp_html}
FASTP_LOG=${fastp_log}
SSC_BAM=${ssc_bam}
SSC_FILTERED_BAM=${ssc_filtered_bam}
SSC_MIN_FAMILY_SIZE=${ssc_min_family_size}
SSC_FAMILY_SIZE_TAG=${ssc_family_size_tag_quoted}

log_step() {
    echo "[$$(date '+%Y-%m-%d %H:%M:%S')] $$*"
}

append_analysis_log() {
    {
        echo ""
        echo "Runtime started: $$(date '+%Y-%m-%d %H:%M:%S')"
        echo "Hostname: $$(hostname)"
        echo "Working directory: $$(pwd)"
        echo "SLURM job ID: $${SLURM_JOB_ID:-not submitted}"
        echo ""
        echo "Runtime program versions:"
        echo "  Java: $$(java -version 2>&1 | head -n 1 || true)"
        echo "  GATK: $$(gatk --version 2>&1 | head -n 1 || true)"
        echo "  BWA: $$(bwa 2>&1 | head -n 3 | tail -n 1 || true)"
        echo "  SAMtools: $$(samtools --version 2>&1 | head -n 1 || true)"
        echo "  BEDTools: $$(bedtools --version 2>&1 | head -n 1 || true)"
        echo "  fastp: $$(fastp --version 2>&1 | head -n 1 || true)"
        echo "  fgbio: $$(java $$JVM_OPTIONS -jar "$$FGBIO_JAR" --version 2>&1 | head -n 1 || true)"
        echo "  pysam: $$(python3 -c 'import pysam; print(pysam.__version__)' 2>&1 || true)"
        echo "  python-Levenshtein: $$(python3 -c 'import Levenshtein; print(Levenshtein.__version__)' 2>&1 || true)"
    } >> "$$ANALYSIS_LOG"
}

finish_analysis_log() {
    local status=$$1
    {
        echo ""
        echo "Runtime finished: $$(date '+%Y-%m-%d %H:%M:%S')"
        echo "Exit status: $$status"
    } >> "$$ANALYSIS_LOG"
}
trap 'finish_analysis_log $$?' EXIT

log_step "Loading Nucleus modules"
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load Java/17.0.15
    module load GATK/4.6.1.0-GCCcore-13.3.0-Java-17
    module load BWA/0.7.19-GCCcore-13.3.0
    module load SAMtools/1.22.1-GCC-13.3.0
    module load BEDTools/2.31.1-GCC-13.3.0
    module load Miniforge3/24.11.3-0
else
    echo "FATAL: module command is not available. Run this job on ERIS Nucleus." >&2
    exit 1
fi

log_step "Initializing conda"
if [ -n "$${EBROOTMINIFORGE3:-}" ] && [ -f "$$EBROOTMINIFORGE3/etc/profile.d/conda.sh" ]; then
    . "$$EBROOTMINIFORGE3/etc/profile.d/conda.sh"
elif [ -n "$${CONDA_PREFIX:-}" ] && [ -f "$$CONDA_PREFIX/etc/profile.d/conda.sh" ]; then
    . "$$CONDA_PREFIX/etc/profile.d/conda.sh"
elif [ -f "$$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    . "$$HOME/miniforge3/etc/profile.d/conda.sh"
else
    echo "FATAL: conda.sh not found after loading Miniforge3." >&2
    exit 1
fi

conda activate fastp-nucleus || {
    echo "FATAL: conda activate fastp-nucleus failed." >&2
    exit 1
}

for tool in gatk bwa samtools bedtools fastp java python3; do
    if ! command -v "$$tool" >/dev/null 2>&1; then
        echo "FATAL: $$tool not found on PATH after module/conda setup." >&2
        exit 1
    fi
done

append_analysis_log

log_step "Preparing directories"
rm -rf "$$TEMP_DIR" "$$TEMP_SSC_DIR" "$$TEMP_DSC_DIR"
mkdir -p "$$TEMP_DIR" "$$TEMP_SSC_DIR" "$$TEMP_DSC_DIR" "$$OUTPUT_SSC_DIR" "$$OUTPUT_DSC_DIR" "$$QC_DIR"

for input_file in "$$R1_FASTQ" "$$R2_FASTQ" "$$REFERENCE_GENOME" "$$FGBIO_JAR" "$$UMI_EXTRACT_SCRIPT" "$$UMI_CATALOG_FILE"; do
    if [ ! -s "$$input_file" ]; then
        echo "FATAL: required input file missing or empty: $$input_file" >&2
        exit 1
    fi
done

log_step "Converting FASTQ to unmapped BAM"
gatk --java-options "$$JVM_OPTIONS" FastqToSam \
    -F1 "$$R1_FASTQ" \
    -F2 "$$R2_FASTQ" \
    -O "$$TEMP_DIR/$${SAMPLE_ID}_unmapped.bam" \
    -SM "$$SAMPLE_ID"

log_step "Classifying and extracting UMIs by species"
python3 "$$UMI_EXTRACT_SCRIPT" \
    --input-bam "$$TEMP_DIR/$${SAMPLE_ID}_unmapped.bam" \
    --output-bam "$$TEMP_DIR/$${SAMPLE_ID}_unmapped_umi_extracted.bam" \
    --work-dir "$$TEMP_DIR" \
    --sample-id "$$SAMPLE_ID" \
    --fgbio-jar "$$FGBIO_JAR" \
    --jvm-options "$$JVM_OPTIONS" \
    --umi-catalog-file "$$UMI_CATALOG_FILE" \
    --three-nt-trim-length "$$THREE_NT_TRIM_LENGTH" \
    | tee -a "$$ANALYSIS_LOG"

log_step "Converting UMI-extracted BAM to FASTQ"
gatk --java-options "$$JVM_OPTIONS" SamToFastq \
    INPUT="$$TEMP_DIR/$${SAMPLE_ID}_unmapped_umi_extracted.bam" \
    FASTQ="$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_R1.fastq" \
    SECOND_END_FASTQ="$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_R2.fastq"

log_step "Running fastp"
fastp \
    -i "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_R1.fastq" \
    -o "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_trimmed_R1.fastq" \
    -I "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_R2.fastq" \
    -O "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_trimmed_R2.fastq" \
    -g -W 5 -q 20 -u 40 -3 -l 75 -c \
    -j "$$FASTP_JSON" \
    -h "$$FASTP_HTML" \
    &> "$$FASTP_LOG"

log_step "Aligning trimmed reads"
bwa mem \
    -t "$$CPU" -K 100000000 \
    -R "@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:$${SAMPLE_ID}\tPU:unit1" \
    -M "$$REFERENCE_GENOME" \
    "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_trimmed_R1.fastq" \
    "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_trimmed_R2.fastq" \
    | samtools view -Sb - > "$$TEMP_DIR/$${SAMPLE_ID}_umi_aligned.bam"

log_step "Merging aligned BAM with UMI metadata"
gatk --java-options "$$JVM_OPTIONS" MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --ALIGNED_BAM "$$TEMP_DIR/$${SAMPLE_ID}_umi_aligned.bam" \
    --UNMAPPED_BAM "$$TEMP_DIR/$${SAMPLE_ID}_unmapped_umi_extracted.bam" \
    --OUTPUT "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_aligned_merged.bam" \
    --REFERENCE_SEQUENCE "$$REFERENCE_GENOME" \
    --SORT_ORDER queryname \
    --ALIGNED_READS_ONLY true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false

log_step "Filtering FR-oriented properly and discordantly paired mapped reads with MAPQ >= 1"
samtools view -h "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_aligned_merged.bam" \
    | awk '/^@/ || (($$2==83 || $$2==99 || $$2==163 || $$2==147 || $$2==81 || $$2==97 || $$2==161 || $$2==145) && $$5>=1)' \
    | samtools view -bh - \
    > "$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_aligned_merged_filtered.bam"

log_step "Grouping SSC reads by UMI"
java $$JVM_OPTIONS -jar "$$FGBIO_JAR" GroupReadsByUmi \
    --input="$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_aligned_merged_filtered.bam" \
    --output="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_grouped.bam" \
    --strategy=adjacency \
    --edits=1 \
    -t RX \
    -f "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_group_data.txt"

log_step "Grouping DSC reads by UMI"
java $$JVM_OPTIONS -jar "$$FGBIO_JAR" GroupReadsByUmi \
    --input="$$TEMP_DIR/$${SAMPLE_ID}_umi_extracted_aligned_merged_filtered.bam" \
    --output="$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_grouped.bam" \
    --strategy=paired \
    --edits=1 \
    -t RX \
    -f "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_group_data.txt"

rm -rf "$$TEMP_DIR"

log_step "Calling SSC consensus reads"
java $$JVM_OPTIONS -jar "$$FGBIO_JAR" CallMolecularConsensusReads \
    --input="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_grouped.bam" \
    --output="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    --error-rate-post-umi 40 \
    --error-rate-pre-umi 45 \
    --min-reads 1 \
    --max-reads 50 \
    --min-input-base-quality 20 \
    --read-name-prefix="consensus"

log_step "Converting SSC consensus BAM to FASTQ"
gatk --java-options "$$JVM_OPTIONS" SamToFastq \
    INPUT="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    FASTQ="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R1.fastq" \
    SECOND_END_FASTQ="$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R2.fastq"

log_step "Realigning SSC consensus reads"
bwa mem \
    -R "@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:$${SAMPLE_ID}\tPU:unit1" \
    -v 3 -Y -M -t "$$CPU" -K 100000000 \
    "$$REFERENCE_GENOME" \
    "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R1.fastq" \
    "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R2.fastq" \
    | samtools view -bh - > "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped_unsort.bam"

log_step "Sorting SSC consensus BAMs by query name"
gatk --java-options "$$JVM_OPTIONS" SortSam \
    -I "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped_unsort.bam" \
    -O "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped.bam" \
    -SO queryname

gatk --java-options "$$JVM_OPTIONS" SortSam \
    -I "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    -O "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_sort.bam" \
    -SO queryname

log_step "Merging SSC consensus alignment"
gatk --java-options "$$JVM_OPTIONS" MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_RETAIN RX \
    --ALIGNED_BAM "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped.bam" \
    --UNMAPPED_BAM "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_sort.bam" \
    --OUTPUT "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_deduped.bam" \
    --REFERENCE_SEQUENCE "$$REFERENCE_GENOME" \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false

log_step "Sorting and indexing SSC BAM"
samtools sort "$$TEMP_SSC_DIR/$${SAMPLE_ID}_umi_deduped.bam" -o "$$SSC_BAM"
samtools index "$$SSC_BAM"

log_step "Generating SSC coverage"
samtools coverage "$$SSC_BAM" \
    | awk 'BEGIN {OFS="\t"} NR==1 {print $$1,$$4,$$6} NR>1 {print $$1,$$4,$$6}' \
    > "$$OUTPUT_SSC_DIR/$${SAMPLE_ID}.SSC_coverage.txt"

log_step "Filtering SSC BAM by family size tag $$SSC_FAMILY_SIZE_TAG > $$SSC_MIN_FAMILY_SIZE"
samtools view -h "$$SSC_BAM" \
    | awk -v min_family_size="$$SSC_MIN_FAMILY_SIZE" -v tag="$$SSC_FAMILY_SIZE_TAG" '
        BEGIN {FS=OFS="\t"}
        /^@/ {print; next}
        {
            keep=0
            for (i=12; i<=NF; i++) {
                split($$i, fields, ":")
                if (fields[1] == tag && fields[2] == "i" && fields[3] > min_family_size) {
                    keep=1
                    break
                }
            }
            if (keep) print
        }
    ' \
    | samtools view -bh - > "$$SSC_FILTERED_BAM"
samtools index "$$SSC_FILTERED_BAM"

log_step "Generating family-size-filtered SSC coverage"
samtools coverage "$$SSC_FILTERED_BAM" \
    | awk 'BEGIN {OFS="\t"} NR==1 {print $$1,$$4,$$6} NR>1 {print $$1,$$4,$$6}' \
    > "$$OUTPUT_SSC_DIR/$${SAMPLE_ID}.SSC_$${SSC_FAMILY_SIZE_TAG}_gt_$${SSC_MIN_FAMILY_SIZE}_coverage.txt"

{
    echo ""
    echo "SSC family-size filter results:"
    echo "  Source BAM: $$SSC_BAM"
    echo "  Filtered BAM: $$SSC_FILTERED_BAM"
    echo "  Filter expression: $$SSC_FAMILY_SIZE_TAG > $$SSC_MIN_FAMILY_SIZE"
    echo "  Filtered read count: $$(samtools view -c "$$SSC_FILTERED_BAM")"
} >> "$$ANALYSIS_LOG"

rm -rf "$$TEMP_SSC_DIR"

log_step "Collecting duplex metrics"
java $$JVM_OPTIONS -jar "$$FGBIO_JAR" CollectDuplexSeqMetrics \
    -i "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_grouped.bam" \
    -o "$$OUTPUT_DSC_DIR/$${SAMPLE_ID}"

log_step "Calling DSC consensus reads"
java $$JVM_OPTIONS -jar "$$FGBIO_JAR" CallDuplexConsensusReads \
    -i "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_grouped.bam" \
    -o "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    --error-rate-post-umi 40 \
    --error-rate-pre-umi 45 \
    --min-reads 1 1 1 \
    --max-reads-per-strand 50 \
    --read-name-prefix="consensus"

log_step "Converting DSC consensus BAM to FASTQ"
gatk --java-options "$$JVM_OPTIONS" SamToFastq \
    INPUT="$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    FASTQ="$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R1.fastq" \
    SECOND_END_FASTQ="$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R2.fastq"

log_step "Realigning DSC consensus reads"
bwa mem \
    -R "@RG\tID:A\tDS:KAPA_TE\tPL:ILLUMINA\tLB:lib1\tSM:$${SAMPLE_ID}\tPU:unit1" \
    -v 3 -Y -M -t "$$CPU" -K 100000000 \
    "$$REFERENCE_GENOME" \
    "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R1.fastq" \
    "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_R2.fastq" \
    | samtools view -bh - > "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped_unsort.bam"

log_step "Sorting DSC consensus BAMs by query name"
gatk --java-options "$$JVM_OPTIONS" SortSam \
    -I "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped_unsort.bam" \
    -O "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped.bam" \
    -SO queryname

gatk --java-options "$$JVM_OPTIONS" SortSam \
    -I "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped.bam" \
    -O "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_sort.bam" \
    -SO queryname

log_step "Merging DSC consensus alignment"
gatk --java-options "$$JVM_OPTIONS" MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_RETAIN RX \
    --ALIGNED_BAM "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_mapped.bam" \
    --UNMAPPED_BAM "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_consensus_unmapped_sort.bam" \
    --OUTPUT "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_deduped.bam" \
    --REFERENCE_SEQUENCE "$$REFERENCE_GENOME" \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false

log_step "Sorting and indexing DSC BAM"
samtools sort "$$TEMP_DSC_DIR/$${SAMPLE_ID}_umi_deduped.bam" -o "$$OUTPUT_DSC_DIR/$${SAMPLE_ID}_umi_dedup_sorted_DSC.bam"
samtools index "$$OUTPUT_DSC_DIR/$${SAMPLE_ID}_umi_dedup_sorted_DSC.bam"

log_step "Generating DSC coverage"
samtools coverage "$$OUTPUT_DSC_DIR/$${SAMPLE_ID}_umi_dedup_sorted_DSC.bam" \
    | awk 'BEGIN {OFS="\t"} NR==1 {print $$1,$$4,$$6} NR>1 {print $$1,$$4,$$6}' \
    > "$$OUTPUT_DSC_DIR/$${SAMPLE_ID}.DSC_coverage.txt"

rm -rf "$$TEMP_DSC_DIR"
log_step "SSC/DSC analysis complete"
''')

    substitutions = {
        "sample_id": args.sample_id,
        "sample_id_quoted": shell_quote(args.sample_id),
        "cpus": str(args.cpus),
        "slurm_mem": args.mem,
        "slurm_time": args.time,
        "partition": args.partition,
        "slurm_stdout": shell_quote(paths["slurm_stdout"]),
        "slurm_stderr": shell_quote(paths["slurm_stderr"]),
        "r1_fastq": shell_quote(paths["r1_fastq"]),
        "r2_fastq": shell_quote(paths["r2_fastq"]),
        "reference_genome": shell_quote(resolve_path(args.reference_genome)),
        "fgbio_jar": shell_quote(resolve_path(args.fgbio_jar)),
        "umi_extract_script": shell_quote(resolve_path(args.umi_extract_script)),
        "umi_catalog_file": shell_quote(resolve_path(args.umi_catalog_file)),
        "three_nt_trim_length": str(args.three_nt_trim_length),
        "jvm_options": shell_quote(args.java_memory),
        "analysis_dir": shell_quote(paths["analysis_dir"]),
        "temp_dir": shell_quote(paths["temp_dir"]),
        "temp_ssc_dir": shell_quote(paths["temp_ssc_dir"]),
        "temp_dsc_dir": shell_quote(paths["temp_dsc_dir"]),
        "output_ssc_dir": shell_quote(paths["output_ssc_dir"]),
        "output_dsc_dir": shell_quote(paths["output_dsc_dir"]),
        "qc_dir": shell_quote(paths["qc_dir"]),
        "analysis_log": shell_quote(paths["analysis_log"]),
        "fastp_json": shell_quote(paths["fastp_json"]),
        "fastp_html": shell_quote(paths["fastp_html"]),
        "fastp_log": shell_quote(paths["fastp_log"]),
        "ssc_bam": shell_quote(paths["ssc_bam"]),
        "ssc_filtered_bam": shell_quote(paths["ssc_filtered_bam"]),
        "ssc_min_family_size": str(args.ssc_min_family_size),
        "ssc_family_size_tag_quoted": shell_quote(args.ssc_family_size_tag),
    }
    return template.substitute(substitutions)


def write_job_script(job_script_path, content):
    job_script_path.write_text(content)
    current_mode = job_script_path.stat().st_mode
    job_script_path.chmod(current_mode | stat.S_IXUSR | stat.S_IXGRP)


def submit_job(job_script_path):
    result = subprocess.run(
        ["sbatch", str(job_script_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.stdout:
        print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())
    if result.returncode != 0:
        raise RuntimeError(f"sbatch failed with exit code {result.returncode}")


def main():
    args = parse_args()
    validate_args(args)
    sample_jobs, skipped_inputs = discover_sample_jobs(args)

    for skipped_input in skipped_inputs:
        print(f"Warning: {skipped_input}")

    if not sample_jobs:
        raise RuntimeError(f"No samples with paired FASTQ inputs found in {build_run_dir(args)}")

    pending_jobs = []
    skipped_complete = []

    for sample_job in sample_jobs:
        validate_sample_inputs(args, sample_job)
        sample_args = build_sample_args(args, sample_job)
        paths = build_paths(sample_args)
        if not args.force and outputs_are_complete(paths):
            skipped_complete.append(sample_job.sample_id)
            continue
        pending_jobs.append((sample_job, sample_args, paths))

    print(f"Discovered {len(sample_jobs)} sample(s) with paired FASTQs in {build_run_dir(args)}")
    if skipped_complete:
        print(f"Skipping {len(skipped_complete)} completed sample(s): {', '.join(skipped_complete)}")
    if not pending_jobs:
        print("No SSC/DSC jobs need to be created or submitted.")
        return

    print(f"Preparing {len(pending_jobs)} SSC/DSC job(s): {', '.join(job[0].sample_id for job in pending_jobs)}")
    if args.dry_run:
        print("Dry run enabled: job scripts and logs will be created, but sbatch submission is skipped.")

    for sample_job, sample_args, paths in pending_jobs:
        create_directories(paths)
        write_analysis_log(sample_args, paths)
        job_script_content = build_job_script(sample_args, paths)
        write_job_script(paths["job_script"], job_script_content)

        print(f"Created job script for {sample_job.sample_id}: {paths['job_script']}")
        print(f"Created analysis log for {sample_job.sample_id}: {paths['analysis_log']}")

        if not args.dry_run:
            submit_job(paths["job_script"])


if __name__ == "__main__":
    main()
