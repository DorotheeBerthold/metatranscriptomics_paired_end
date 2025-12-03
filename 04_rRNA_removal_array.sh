#!/bin/bash
#SBATCH --job-name=rRNA_fil
#SBATCH --output=rRNA_%A_%array.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=50G
#SBATCH --array=1-16       # adjust to total number of samples
set -euo pipefail

# Load modules
module load vsearch
module load python

# Infernal path
export PATH="$HOME/infernal-1.1.5/bin:$PATH"

# Directories
QUAL_UNIVEC_DIR="quality_vector_rm_reads"
RRNA_DIR="rRNA_rm_OMM12_reads"
mkdir -p "$RRNA_DIR"

# Rfam CM DB
RFAM_DB="OMM12_rRNA_DB/OMM12_rRNA.cm"

# Sample list
SAMPLE_LIST="samples.txt"   

# Pick sample for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
fastq="${QUAL_UNIVEC_DIR}/${sample}"

# --- Remove .fastq extension for output prefix ---
prefix="${sample%.fastq}"

# --- Output files ---
fasta="${RRNA_DIR}/${prefix}.fasta"
infernal_tbl="${RRNA_DIR}/${prefix}_rRNA.infernalout"
infernal_log="${RRNA_DIR}/${prefix}_rRNA.log"

echo "Processing sample: $sample"

# 1. Convert FASTQ -> FASTA
vsearch --fastq_filter "$fastq" --fastaout "$fasta"

# 2. Run Infernal cmsearch
cmsearch --notextw \
        --cpu 8 \
         -o "$infernal_log" \
         --tblout "$infernal_tbl" \
         --anytrunc --rfam -E 0.001 \
         "$RFAM_DB" "$fasta"

echo "Finished sample: $sample"