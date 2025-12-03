#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc_%j.log
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --cpus-per-task=4
#SBATCH --array=0-15

set -euo pipefail


# Load FastQC module (check available versions with 'module avail fastqc')
module load fastqc

# Directories
RAW_DIR="trimmed_reads"
FASTQC_DIR="trimmed_fastqc"

mkdir -p "$FASTQC_DIR"

# Get sample from list
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" samples_trimmed_fastqc.txt)
input_file="${RAW_DIR}/${sample}"

# Run FastQC
fastqc "$input_file" -o "$FASTQC_DIR"