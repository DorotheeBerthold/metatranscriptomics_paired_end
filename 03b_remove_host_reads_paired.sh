#!/bin/bash
#SBATCH --job-name=remove_host_counted
#SBATCH --array=0-5                              
#SBATCH --cpus-per-task=4                        
#SBATCH --mem=10G                                
#SBATCH --output=logs/host_rm_%A_%a.out
#SBATCH --error=logs/host_rm_%A_%a.err

# Directories
RAW_DIR="quality_vector_rm_reads"
OUTPUT_DIR="clean_reads"
STATS_DIR="stats"
DB="mouse_cds.fa"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$STATS_DIR"
mkdir -p logs

# Initialize the CSV header if it doesn't exist
SUMMARY_FILE="${STATS_DIR}/cleaning_summary.csv"
if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Sample,Raw_Pairs,BWA_Clean_Pairs,BWA_Dropped,Final_Pairs,BLAT_Dropped,Total_Dropped" > "$SUMMARY_FILE"
fi

# Helper function to count reads (lines divided by 4)
count_reads() {
    local file="$1"
    echo $(($(wc -l < "$file") / 4))
}

# 1. Get the sample prefix
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" samples.txt)
fwd="${RAW_DIR}/${sample}_R1_univec_clean.fastq"
rev="${RAW_DIR}/${sample}_R2_univec_clean.fastq"
sam_out="${OUTPUT_DIR}/${sample}_mouse_bwa.sam"
bam_out="${OUTPUT_DIR}/${sample}_mouse_bwa.bam"

echo "Processing sample: ${sample}"

# Count initial raw pairs (Just need to count R1 since they are paired)
raw_pairs=$(count_reads "$fwd")

# 2. BWA Alignment & Filtering
bwa mem -t $SLURM_CPUS_PER_TASK "$DB" "$fwd" "$rev" > "$sam_out"
samtools view -bS "$sam_out" > "$bam_out"

bwa_clean_r1="${OUTPUT_DIR}/${sample}_R1_bwa_clean.fastq"
bwa_clean_r2="${OUTPUT_DIR}/${sample}_R2_bwa_clean.fastq"

samtools fastq -f 12 -1 "$bwa_clean_r1" -2 "$bwa_clean_r2" -0 /dev/null -s /dev/null "$bam_out"
rm "$sam_out" "$bam_out"

# Count pairs after BWA
bwa_pairs=$(count_reads "$bwa_clean_r1")
bwa_dropped=$((raw_pairs - bwa_pairs))

# ===========================================================================================
# 4. Secondary Alignment with BLAT (Paired-End Aware)
# ===========================================================================================

final_r1="${OUTPUT_DIR}/${sample}_R1_final.fastq"
final_r2="${OUTPUT_DIR}/${sample}_R2_final.fastq"

# Generate FASTA files for BLAT
vsearch --fastq_filter "$bwa_clean_r1" --fastaout "${OUTPUT_DIR}/${sample}_R1.fasta"
vsearch --fastq_filter "$bwa_clean_r2" --fastaout "${OUTPUT_DIR}/${sample}_R2.fasta"

# Run BLAT on R1 and R2
blat -noHead -minIdentity=90 -minScore=65 "$DB" "${OUTPUT_DIR}/${sample}_R1.fasta" -fine \
     -q=rna -t=dna -out=blast8 "${OUTPUT_DIR}/${sample}_R1.blatout"

blat -noHead -minIdentity=90 -minScore=65 "$DB" "${OUTPUT_DIR}/${sample}_R2.fasta" -fine \
     -q=rna -t=dna -out=blast8 "${OUTPUT_DIR}/${sample}_R2.blatout"

# Run the paired Python script to filter them synchronously
./1_BLAT_Filter_DB_paired.py \
    "$bwa_clean_r1" "$bwa_clean_r2" \
    "${OUTPUT_DIR}/${sample}_R1.blatout" "${OUTPUT_DIR}/${sample}_R2.blatout" \
    "$final_r1" "$final_r2" \
    "${OUTPUT_DIR}/${sample}_R1_blat_contam.fastq" "${OUTPUT_DIR}/${sample}_R2_blat_contam.fastq"

# Cleanup intermediate files
rm "${OUTPUT_DIR}/${sample}_R1.fasta" "${OUTPUT_DIR}/${sample}_R2.fasta"
rm "${OUTPUT_DIR}/${sample}_R1.blatout" "${OUTPUT_DIR}/${sample}_R2.blatout"
rm "$bwa_clean_r1" "$bwa_clean_r2"