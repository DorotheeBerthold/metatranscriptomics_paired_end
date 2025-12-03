#!/bin/bash
#SBATCH --job-name=trim_vsearch_qc
#SBATCH --output=trim_vsearch_qc_%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=16000
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4

set -euo pipefail

# Load required modules
module load bwa
module load samtools
module load vsearch

# Define paths
QUAL_DIR="trimmed_reads"
UNIVEC_DIR="UniVec_Core_db"
QUAL_UNIVEC_DIR="quality_vector_rm_reads"

mkdir -p "$QUAL_UNIVEC_DIR"

# make sure samples.txt only contains names without any R1/R2 etc

# Get base sample name from array file
base=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples2.txt)

# Construct full paths
fwd="${QUAL_DIR}/${base}_R1_paired.fastq"
rev="${QUAL_DIR}/${base}_R2_paired.fastq"


# Output prefix
out="${QUAL_UNIVEC_DIR}/${base}_univec"

# Align reads to UniVec
bwa mem -t 4 "${UNIVEC_DIR}/UniVec_Core" "$fwd" "$rev" > "${out}.sam"

# Convert SAM to BAM
samtools view -bS "${out}.sam" > "${out}.bam"

# Extract reads that did NOT map to UniVec (clean reads)
samtools fastq -n -f 12 \
    -1 "${QUAL_UNIVEC_DIR}/${base}_R1_univec_clean.fastq" \
    -2 "${QUAL_UNIVEC_DIR}/${base}_R2_univec_clean.fastq" \
    "${out}.bam"

# Optionally extract reads that DID map (contaminants)
samtools fastq -n -f 3 \
    -1 "${QUAL_UNIVEC_DIR}/${base}_R1_univec_contaminants.fastq" \
    -2 "${QUAL_UNIVEC_DIR}/${base}_R2_univec_contaminants.fastq" \
    "${out}.bam"

# Clean up
rm "${out}.sam" "${out}.bam"

