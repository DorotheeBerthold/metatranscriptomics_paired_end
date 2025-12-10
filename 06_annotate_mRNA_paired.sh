#!/bin/bash
#SBATCH --job-name=annotate_mRNA
#SBATCH --output=annotate_mRNA_%A_%a.log
#SBATCH --array=0-7      # Adjust after checking number of files
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16G
set -euo pipefail

module load bwa
module load samtools
module load vsearch

############### RUN BEFORE SUBMITTING TO SBATCH #############################
#bwa index -a bwtsw "$REF"
#############################################################################

# -------- INPUTS --------
READ_DIR="rRNA_filtered"
REF="OMM12_pseudogenome.fasta"        # Your pseudogenome

# Get list of R1 files only
FILES=($(ls ${READ_DIR}/*_R1_unique_mRNA.fastq))
FILE_R1=${FILES[$SLURM_ARRAY_TASK_ID]}

# Infer R2 file
FILE_R2=${FILE_R1/_R1_/_R2_}

# Extract sample name
BASENAME=$(basename "$FILE_R1")
SAMPLE=${BASENAME%_R1_unique_mRNA.fastq}


# Index the microbial genome database
samtools faidx "$REF"

echo "Processing sample: $SAMPLE"

# -------- OUTPUT DIRS --------
OUTDIR="annotation_results/${SAMPLE}"
mkdir -p ${OUTDIR}


# -------- STEP 2 — Annotate reads --------
echo "Running BWA MEM on reads..."
bwa mem -t 4 "$REF" "$FILE_R1" "$FILE_R2" \
    > "${OUTDIR}/${SAMPLE}_annotation_bwa.sam"


echo "Done with $SAMPLE"
