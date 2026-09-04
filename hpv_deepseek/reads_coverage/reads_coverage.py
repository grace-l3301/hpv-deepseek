#!/usr/bin/env python3
"""
reads_coverage.py

Generate and submit SLURM job scripts for UMI deduplication and coverage analysis pipeline.
Dynamically detects samples from directory structure instead of using fixed count.

Author: Samuli Eldfors
"""

__version__ = "9.1.0"

import os
import subprocess
import argparse
import glob
import re
from datetime import datetime

# Define variables
BASE_DIR = "/data/fadenlab/processed"
RUN_DATE = "10_27_2025-test"

REFERENCE_GENOME = "/data/fadenlab/reference/HPV-DeepSeek/HPV_hg38/HPV_hg38.fasta"
FGBIO_JAR = "/data/fadenlab/bin/fgbio/fgbio-2.1.0.jar"
THREAD_COUNT = 2
JOB_SCRIPT_DIR = os.path.join(BASE_DIR, "UMI", RUN_DATE, "job_scripts")
LOG_DIR = os.path.join(BASE_DIR, "UMI", RUN_DATE)
LOG_FILE = os.path.join(LOG_DIR, "pipeline_log.txt")
JVM_OPTIONS = "-Xmx32g"

# Create job script directory and log directory if they don't exist
os.makedirs(JOB_SCRIPT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
print(f"Job script directory created: {JOB_SCRIPT_DIR}")
print(f"Log directory created: {LOG_DIR}")

# Parse arguments
parser = argparse.ArgumentParser(description="Generate and (optionally) submit SLURM job scripts for UMI pipeline (HPV_EXP01).")
parser.add_argument("--dry-run", action="store_true", help="Create job scripts but do not submit them with sbatch.")
args = parser.parse_args()
DRY_RUN = args.dry_run
if DRY_RUN:
    print("Dry run enabled: job scripts will be created but not submitted (sbatch skipped).")

# Auto-detect samples from directory structure
# Look for directories matching pattern: {lib_id}_S{rep_id}_L001
sample_dirs = glob.glob(os.path.join(BASE_DIR, "UMI", RUN_DATE, "*_S*_L001"))
sample_pattern = re.compile(r'(\d+_S\d+_L001)$')

samples = []
for dir_path in sample_dirs:
    match = sample_pattern.search(dir_path)
    if match:
        sample_name = match.group(1)
        # Verify input directory exists with files
        input_dir = os.path.join(dir_path, "input")
        if os.path.isdir(input_dir):
            r1_file = os.path.join(input_dir, f"{sample_name}_R1_001.fastq.gz")
            r2_file = os.path.join(input_dir, f"{sample_name}_R2_001.fastq.gz")
            if os.path.exists(r1_file) and os.path.exists(r2_file):
                samples.append(sample_name)
            else:
                print(f"Warning: Input files missing for {sample_name}")

samples.sort()
print(f"Found {len(samples)} samples to process")
if len(samples) == 0:
    print("Error: No samples found. Exiting.")
    exit(1)

# Write log file with program versions and other details
with open(LOG_FILE, "w") as log_file:
    log_file.write(f"Script name: reads_coverage-v{__version__}.py\n")
    log_file.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log_file.write(f"Base directory: {BASE_DIR}\n")
    log_file.write(f"Run date: {RUN_DATE}\n")
    log_file.write(f"Reference genome: {REFERENCE_GENOME}\n")
    log_file.write(f"Number of samples: {len(samples)}\n")
    log_file.write(f"Thread count: {THREAD_COUNT}\n")
    log_file.write("Cluster: ERIS Nucleus\n")
    log_file.write("Program versions:\n")
    log_file.write("  - GATK: 4.6.1.0\n")
    log_file.write("  - BWA: 0.7.19\n")
    log_file.write("  - SAMtools: 1.22.1\n")
    log_file.write("  - BEDTools: 2.31.1\n")
    log_file.write("  - FGBIO: 2.1.0\n")
    log_file.write("  - Fastp: (version not specified, but activated via conda environment `fastp`)\n")
    log_file.write(f"\nSamples to process:\n")
    for sample in samples:
        log_file.write(f"  - {sample}\n")
    log_file.write("\n")

# Template for the SLURM job script
job_script_template = """#!/bin/bash
#SBATCH --job-name=umi_pipeline_{sample_name}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=13
#SBATCH --mem=62G
#SBATCH --time=24:00:00
#SBATCH --output={base_dir}/UMI/{run_date}/{sample_name}/output/umi_pipeline_%j.out
#SBATCH --error={base_dir}/UMI/{run_date}/{sample_name}/output/umi_pipeline_%j.err

# Abort on any error, undefined variable, or pipe failure
set -euo pipefail

# ── Module initialization ──────────────────────────────────────────────
# Purge inherited module state from login environment to ensure a clean
# and reproducible software stack on Nucleus compute nodes.
module purge
module load GATK/4.6.1.0-GCCcore-13.3.0-Java-17
module load BWA/0.7.19-GCCcore-13.3.0
module load SAMtools/1.22.1-GCC-13.3.0
module load BEDTools/2.31.1-GCC-13.3.0
module load Miniforge3/24.11.3-0

# ── Conda initialization ──────────────────────────────────────────────
# In a non-interactive Slurm batch shell, conda's shell functions are
# not available by default.  We must explicitly source conda.sh from the
# Miniforge3 installation to define `conda activate` as a shell function.
# On HPC we do NOT rely on `conda init` (which modifies ~/.bashrc) because
# that approach is fragile across different clusters and login nodes.

# EBROOTMINIFORGE3 is set by the EasyBuild Miniforge3 module.
if [ -n "${{EBROOTMINIFORGE3:-}}" ] && [ -f "$EBROOTMINIFORGE3/etc/profile.d/conda.sh" ]; then
    . "$EBROOTMINIFORGE3/etc/profile.d/conda.sh"
else
    echo "FATAL: EBROOTMINIFORGE3 is not set or conda.sh not found after 'module load Miniforge3'." >&2
    echo "       Verify Miniforge3 module with: module avail Miniforge3" >&2
    exit 1
fi

# Diagnostics – logged so we can verify correct conda in job output
echo "--- conda diagnostics ---"
which conda
type conda
echo "CONDA_PREFIX=${{CONDA_PREFIX:-unset}}"

# ── Activate fastp environment ────────────────────────────────────────
conda activate fastp-nucleus
echo "Activated conda env: CONDA_PREFIX=$CONDA_PREFIX"

# Verify fastp is actually available after activation
if ! command -v fastp >/dev/null 2>&1; then
    echo "FATAL: fastp not found on PATH after 'conda activate fastp-nucleus'." >&2
    echo "       The environment may need to be rebuilt on Nucleus." >&2
    echo "       Try: conda env create -f /data/fadenlab/scripts/scripts-nucleus/envs/fastp.yml" >&2
    exit 1
fi
fastp --version 2>&1 || true
echo "--- end conda diagnostics ---"

# Job information
echo "Running job for sample {sample_name}"
echo "Base directory: {base_dir}"
echo "Run date: {run_date}"
echo "Reference genome: {reference_genome}"

# Create temp directory
mkdir -p {base_dir}/UMI/{run_date}/{sample_name}/temp

# Check if input files exist
if [ ! -f {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R1_001.fastq.gz ]; then
    echo "Error: Input file {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R1_001.fastq.gz not found."
    exit 1
fi
if [ ! -f {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R2_001.fastq.gz ]; then
    echo "Error: Input file {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R2_001.fastq.gz not found."
    exit 1
fi

# Convert FastQ to BAM
gatk --java-options "{jvm_options}" FastqToSam \\
-F1 {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R1_001.fastq.gz \\
-F2 {base_dir}/UMI/{run_date}/{sample_name}/input/{sample_name}_R2_001.fastq.gz \\
-O {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped.bam \\
-SM {sample_name}

# Check if the output file exists and is non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped.bam ]; then
    echo "Error: {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped.bam not found or empty."
    exit 1
fi

# Extract UMIs from BAM 
java {jvm_options} -jar {fgbio_jar} ExtractUmisFromBam \\
-i {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped.bam \\
-o {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped_umi_extracted.bam \\
-r 3M3S+T 3M3S+T \\
-t RX \\
-a true

# Check if the output file exists and is non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped_umi_extracted.bam ]; then
    echo "Error: {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped_umi_extracted.bam not found or empty."
    exit 1
fi

# Adapter trimming and quality filtering; first need to convert BAM back to Fastq
gatk --java-options "{jvm_options}" SamToFastq \\
INPUT={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped_umi_extracted.bam \\
FASTQ={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R1.fastq \\
SECOND_END_FASTQ={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R2.fastq

# Check if the output files exist and are non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R1.fastq ] || [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R2.fastq ]; then
    echo "Error: Fastq files not found or empty."
    exit 1
fi

# Adapter trimming and quality filtering
fastp \\
-i {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R1.fastq \\
-o {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R1.fastq \\
-I {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_R2.fastq \\
-O {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R2.fastq \\
-g -W 5 -q 20 -u 40 -3 -l 75 -c \\
-j {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_fastp.json \\
-h {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_fastp.html \\
&> {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_fastp.log

# Check if the output files exist and are non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R1.fastq ] || [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R2.fastq ]; then
    echo "Error: Trimmed Fastq files not found or empty."
    exit 1
fi

# Map reads to HPV_hg38 genome
# Write to SAM first to avoid pipe issues, then convert to BAM
bwa mem \\
-t {thread_count} -K 100000000 \\
-R '@RG\\tID:A\\tDS:KAPA_TE\\tPL:ILLUMINA\\tLB:lib1\\tSM:sample1\\tPU:unit1' \\
-M {reference_genome} \\
{base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R1.fastq \\
{base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_trimmed_R2.fastq \\
-o {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.sam

# Convert SAM to BAM
samtools view -Sb {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.sam > {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.bam
rm {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.sam

# Check if the output file exists and is non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.bam ]; then
    echo "Error: {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.bam not found or empty."
    exit 1
fi

# Add UMI info (extracted from before) to the aligned BAMs
gatk MergeBamAlignment \\
--ATTRIBUTES_TO_RETAIN X0 \\
--ATTRIBUTES_TO_REMOVE NM \\
--ATTRIBUTES_TO_REMOVE MD \\
--ALIGNED_BAM {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_aligned.bam \\
--UNMAPPED_BAM {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_unmapped_umi_extracted.bam \\
--OUTPUT {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged.bam \\
--REFERENCE_SEQUENCE {reference_genome} \\
--SORT_ORDER queryname \\
--ALIGNED_READS_ONLY true \\
--MAX_INSERTIONS_OR_DELETIONS -1 \\
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
--ALIGNER_PROPER_PAIR_FLAGS true \\
--CLIP_OVERLAPPING_READS false

# Check if the output file exists and is non-empty
if [ ! -s {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged.bam ]; then
    echo "Error: {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged.bam not found or empty."
    exit 1
fi

# Filter reads with samtools view
samtools view -f 2 -q 1 -bh {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged.bam > {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged_filtered.bam

# Group reads from the same UMI
java {jvm_options} -jar {fgbio_jar} GroupReadsByUmi \\
--input={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_extracted_aligned_merged_filtered.bam \\
--output={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_grouped.bam \\
--strategy=adjacency \\
--edits=1 \\
-t RX \\
-f {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_group_data.txt

# Generate consensus sequences from groups of reads that all come from the same UMI
java {jvm_options} -jar {fgbio_jar} CallMolecularConsensusReads \\
--input={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_grouped.bam \\
--output={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped.bam \\
--error-rate-post-umi 40 \\
--error-rate-pre-umi 45 \\
--output-per-base-tags false \\
--min-reads 1 \\
--max-reads 50 \\
--min-input-base-quality 20 \\
--read-name-prefix='consensus'

# Convert consensus BAM back to fastq so that we can reobtain alignment coordinate information
gatk SamToFastq \\
INPUT={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped.bam \\
FASTQ={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_R1.fastq \\
SECOND_END_FASTQ={base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_R2.fastq

# Realign consensus reads to reference genome
# Write to SAM first to avoid pipe issues, then convert to BAM
bwa mem \\
-t {thread_count} \\
-R '@RG\\tID:A\\tDS:KAPA_TE\\tPL:ILLUMINA\\tLB:lib1\\tSM:sample1\\tPU:unit1' \\
-v 3 -Y -M -K 100000000 \\
{reference_genome} \\
{base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_R1.fastq \\
{base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_R2.fastq \\
-o {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped_unsort.sam

# Convert SAM to BAM
samtools view -bh {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped_unsort.sam > {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped_unsort.bam
rm {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped_unsort.sam

# Sort both unmapped and mapped consensus bams by queryname
gatk --java-options "{jvm_options}" SortSam \\
-I {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped_unsort.bam \\
-O {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped.bam \\
-SO queryname

gatk --java-options "{jvm_options}" SortSam \\
-I {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped.bam \\
-O {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_sort.bam \\
-SO queryname

# Add UMI info (from before) to consensus reads in BAM
gatk --java-options "{jvm_options}" MergeBamAlignment \\
--ATTRIBUTES_TO_RETAIN X0 \\
--ATTRIBUTES_TO_RETAIN RX \\
--ALIGNED_BAM {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_mapped.bam \\
--UNMAPPED_BAM {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_consensus_unmapped_sort.bam \\
--OUTPUT {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_deduped.bam \\
--REFERENCE_SEQUENCE {reference_genome} \\
--SORT_ORDER coordinate \\
--ADD_MATE_CIGAR true \\
--MAX_INSERTIONS_OR_DELETIONS -1 \\
--PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
--ALIGNER_PROPER_PAIR_FLAGS true \\
--CLIP_OVERLAPPING_READS false

# Sort and index deduplicated BAM
samtools sort {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_umi_deduped.bam -o {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_umi_dedup_sorted.bam
samtools index {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_umi_dedup_sorted.bam

# Generate coverage file for initial results
samtools coverage {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_umi_dedup_sorted.bam | awk '{{{{if(NR==1){{{{printf "%s\\t%s\\t%s\\n",$1,$4,$6}}}} else {{{{printf "%s\\t%s\\t%s\\n",$1,$4,$6}}}}}}}}' > {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}.coverage.txt

# Filter reads that map to HPV31_Ref:4093-4137
bedtools intersect -abam {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}_umi_dedup_sorted.bam \\
-b <(echo -e "HPV31_Ref\\t4093\\t4137") -v \\
> {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_filtered.bam
samtools index {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_filtered.bam

# Generate coverage file for filtered results
samtools coverage {base_dir}/UMI/{run_date}/{sample_name}/temp/{sample_name}_filtered.bam | \\
awk '{{{{if(NR==1){{{{printf "%s\\t%s\\t%s\\n",$1,$4,$6}}}} else {{{{printf "%s\\t%s\\t%s\\n",$1,$4,$6}}}}}}}}' > {base_dir}/UMI/{run_date}/{sample_name}/output/{sample_name}.filtered.coverage.txt

# Remove temp directory
rm -r {base_dir}/UMI/{run_date}/{sample_name}/temp

"""

# Loop through each detected sample and create job scripts
for sample_name in samples:
    print(f"Processing sample: {sample_name}")
    
    job_script_content = job_script_template.format(
        sample_name=sample_name,
        base_dir=BASE_DIR,
        run_date=RUN_DATE,
        reference_genome=REFERENCE_GENOME,
        fgbio_jar=FGBIO_JAR,
        thread_count=THREAD_COUNT,
        jvm_options=JVM_OPTIONS
    )
    
    job_script_path = os.path.join(JOB_SCRIPT_DIR, f"job_{sample_name}.sh")
    
    # Write the job script to a file
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    
    # Print the path of the job script for debugging
    print(f"Created job script: {job_script_path}")
    
    # Submit the job script to the cluster unless dry-run is enabled
    if DRY_RUN:
        print(f"Dry run: would submit job script with sbatch: {job_script_path}")
    else:
        from subprocess import PIPE
        result = subprocess.run(["sbatch", job_script_path], stdout=PIPE, stderr=PIPE, universal_newlines=True)
        print("Submitted job for {}: {}".format(sample_name, result.stdout.strip()))
        if result.stderr:
            print("Error submitting job for {}: {}".format(sample_name, result.stderr.strip()))

print(f"\nJob scripts created for {len(samples)} samples.")
if not DRY_RUN:
    print("All jobs submitted successfully.")
