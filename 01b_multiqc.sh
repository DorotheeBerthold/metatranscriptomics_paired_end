#!/bin/bash
#SBATCH --job-name=multiqc
#SBATCH --output=multiqc_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2

# Activate multiqc env
source ~/multiqc/bin/activate

FASTQC_DIR="trimmed_fastqc"

multiqc "$FASTQC_DIR" -o "$FASTQC_DIR"/summary_multiqc

