#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=count_%A_%a.log
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-7   # <-- adjust automatically based on number of sample folders

# LOAD MODULES
module load subread
module load samtools

# INPUT LISTx
SAMPLE_DIRS=($(ls -d annotation_results/*))

# SELECT CURRENT SAMPLE
SAMPLE=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}
echo "Processing sample directory: $SAMPLE"

# IDENTIFY INPUT .sam FILE
SAM=$(ls $SAMPLE/*.sam)
echo "Found SAM file: $SAM"

# OUTPUT PATHS
BAM=${SAM%.sam}.bam
SORTED_BAM=${SAM%.sam}_sorted.bam
COUNTS=${SAMPLE}/gene_counts.txt

# REFERENCE ANNOTATION
ANNOTATION=OMM12_pseudogenome_genes2.gff3

###############
# STEP 1: SAM → BAM + SORT
###############
echo "Converting SAM to BAM..."
samtools view -@ 4 -bS "$SAM" > "$BAM"

echo "Sorting BAM..."
samtools sort -@ 4 -o "$SORTED_BAM" "$BAM"

samtools index "$SORTED_BAM"

###############
# STEP 2: HTSEQ COUNT
###############
echo "Running featureCounts..."
featureCounts \
    -T 4 \
    -p \
    -t gene \
    -g locus_tag \
    -a "$ANNOTATION" \
    -o "$COUNTS" \
    "$SORTED_BAM"

echo "DONE for $SAMPLE"