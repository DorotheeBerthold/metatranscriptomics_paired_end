#!/bin/bash
#SBATCH --job-name=trim_vsearch_qc
#SBATCH --output=trim_vsearch_qc_%j.log
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4

set -euo pipefail

module load trimmomatic
module load vsearch

### Define paths
RAW_DIR="raw_files_combined"
TRIM_DIR="trimmed_reads"
ADAPTERS="adapter/NexteraPE-PE.fa"


mkdir -p "$TRIM_DIR" "$QUAL_DIR"

sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" samples.txt)

fwd="${RAW_DIR}/${sample}_R1_001.fastq"
rev="${RAW_DIR}/${sample}_R2_001.fastq"

# run Trimmomatic ONCE per array job
trimmomatic PE -threads 8 \
    "$fwd" "$rev" \
    "${TRIM_DIR}/${sample}_R1_paired.fastq" "${TRIM_DIR}/${sample}_R1_unpaired.fastq" \
    "${TRIM_DIR}/${sample}_R2_paired.fastq" "${TRIM_DIR}/${sample}_R2_unpaired.fastq" \
     HEADCROP:10 \
     ILLUMINACLIP:${ADAPTERS}:2:30:10 \
     LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50


echo "=== Done with sample ${sample} ==="