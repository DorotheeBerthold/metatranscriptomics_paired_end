#!/bin/bash
#SBATCH --job-name=rRNA_filter
#SBATCH --output=rRNAfilter_%A_%a.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-8
set -euo pipefail

module load python/3.12.8


FILTER_SCRIPT="2_Infernal_Filter_paired.py" 
chmod +x "$FILTER_SCRIPT"

QUAL_UNIVEC_DIR="quality_vector_rm_reads"
RRNA_DIR="rRNA_rm_OMM12_reads"
OUTDIR="rRNA_filtered"
mkdir -p "$OUTDIR"

# Sample list without any extension (R1/R2)
SAMPLE_LIST="samples2.txt"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")

# Input FASTQs
R1="${QUAL_UNIVEC_DIR}/${sample}_R1_univec_clean.fastq"
R2="${QUAL_UNIVEC_DIR}/${sample}_R2_univec_clean.fastq"

# Infernal output (already computed)
INF_R1="${RRNA_DIR}/${sample}_R1_univec_clean_rRNA.infernalout"
INF_R2="${RRNA_DIR}/${sample}_R2_univec_clean_rRNA.infernalout"

# Combined rRNA ID list
COMBINED_IDS="${OUTDIR}/${sample}_rRNA_ids.txt"

# Output FASTQs
mRNA_R1="${OUTDIR}/${sample}_R1_unique_mRNA.fastq"
mRNA_R2="${OUTDIR}/${sample}_R2_unique_mRNA.fastq"
rRNA_R1="${OUTDIR}/${sample}_R1_unique_rRNA.fastq"
rRNA_R2="${OUTDIR}/${sample}_R2_unique_rRNA.fastq"

# Skip existing
if [[ -f "$mRNA_R1" && -f "$mRNA_R2" ]]; then
    echo "Skipping $sample — already done"
    exit 0
fi

echo "Processing $sample"

##############################################
# Step 1 — Extract rRNA IDs from Infernal
##############################################
grep -v '^#' "$INF_R1" | awk '{print $1}' > "$COMBINED_IDS"
grep -v '^#' "$INF_R2" | awk '{print $1}' >> "$COMBINED_IDS"

sort -u "$COMBINED_IDS" -o "$COMBINED_IDS"

##############################################
# Step 2 — Run paired-end rRNA filtering
##############################################
./"$FILTER_SCRIPT" \
    "$R1" "$R2" \
    "$COMBINED_IDS" \
    "$mRNA_R1" "$mRNA_R2" \
    "$rRNA_R1" "$rRNA_R2"

echo "Finished $sample"
