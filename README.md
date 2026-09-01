# Metatranscriptomics Paired-End Analysis Pipeline

A comprehensive workflow for processing and analyzing paired-end metatranscriptomic data. This pipeline performs quality control, read preprocessing, contamination removal, taxonomic/functional annotation, and statistical analysis to identify differentially expressed genes across conditions. It has been optimised for paired-end data processed on a high-performance cluster (HPC) to submit via slurm.

---

## Summary of Main Scripts

### Pre-processing & Quality Control
- **01_fastqc.sh** — Run FastQC quality control on trimmed reads
- **01b_multiqc.sh** — Aggregate FastQC results with MultiQC
- **02_processing_reads_paired.sh** — Trim adapters and quality-filter reads using Trimmomatic
- **03_remove_vector_contamination_paired.sh** — Remove UniVec database contamination using BWA/VSEARCH
- **03b_remove_host_reads_paired.sh** — Remove host/non-target reads

### rRNA & Annotation Database Preparation
- **04_rRNA_removal_array.sh** — Filter ribosomal RNA using Infernal and custom rRNA database
- **05_rRNA_filter_paired.sh** — Create final filtered mRNA dataset
- **05b_create_rRNA_overview_graphs.R** — Visualize rRNA removal statistics
- **rRNA_DB/build_rRNA_DB.sh** — Build custom Infernal rRNA database
- **pseudogenome_creation_circular.py** — Create pseudogenome from microbial genome sequences (GBFF format)
- **1_BLAT_Filter_DB_paired.py** — Filter contaminated reads using BLAT results
- **2_Infernal_Filter_paired.py** — Filter contaminated reads using Infernal results

### Read Annotation & Quantification
- **06_annotate_mRNA_paired.sh** — Map mRNA reads to pseudogenome using BWA
- **06b_annotate_rRNA_paired.sh** — Map rRNA reads to pseudogenome
- **07_count_reads_paired.sh** — Count reads mapped to each gene using featureCounts

### Statistical Analysis & Visualization
- **08_merge_counts.R** — Merge read count data from all samples
- **09_PCA.R** — Perform PCA analysis and visualization
- **10_DESeq2_taxon_normalised.R** — Differential expression analysis using DESeq2 (taxon-normalized)
- **11a_functional_annotation_genomes.R** — Annotate genomes with functional information
  - **run_cayman.sh** — Annotate protein sequences with CAYMAN (functional annotation)
  - **eggnog_map_loop.sh** - Annotate protein sequences with eggnog
- **11b_merge_functional_annotations_DESEq2.R** — Merge functional annotations with DESeq2 results
- **11c_extract_cazyme_protein_sequences.R** — Extract CAZyme sequences
- **12_volcanoplots.R** — Generate volcano plots for DESeq2 results
- **13a_KEGG_pathway_enrichment.R** — KEGG pathway enrichment analysis
- **13b_KEGG_pathview_maps.R** — Create KEGG pathway visualization maps
- **13c_GSEA_enrichment.R** — Gene set enrichment analysis

---

## Step-by-Step Pipeline Guide

### **PHASE 0: Setup & Database Preparation**

#### 0.1 Prepare Directory Structure
```bash
# Create main directories for the pipeline
mkdir -p raw_files_combined
mkdir -p trimmed_reads
mkdir -p trimmed_fastqc
mkdir -p quality_vector_rm_reads
mkdir -p rRNA_rm_OMM12_reads
mkdir -p rRNA_filtered
mkdir -p annotation_results
mkdir -p results
mkdir -p OMM_gbff  # For GenBank files
```

#### 0.2 Create Sample File
Create `samples.txt` with one sample name per line (no R1/R2 suffixes):
```
Sample_S01
Sample_S02
Sample_S03
...
```

#### 0.3 Prepare Adapter File
Download or create an adapter FASTA file (e.g., `adapter/NexteraPE-PE.fa` for Illumina NextEra).

#### 0.4 Build rRNA Database (if needed)
```bash
cd rRNA_DB
bash build_rRNA_DB.sh
# This creates OMM12_rRNA.cm database using Infernal
cd ..
```

#### 0.5 Create Pseudogenome from Microbial Genomes
1. Place GenBank files (.gbff) in `OMM_gbff/` directory
2. Run the pseudogenome creation script:
```bash
python pseudogenome_creation_circular.py
# Generates: OMM12_pseudogenome.fasta, OMM12_pseudogenome.gff3
```

#### 0.6 Index the Pseudogenome
```bash
bwa index -a bwtsw OMM12_pseudogenome.fasta
samtools faidx OMM12_pseudogenome.fasta
```

---

### **PHASE 1: Quality Control & Trimming**

#### 1.1 Run Initial FastQC (Optional)
```bash
# Create sample list for raw reads
ls raw_files_combined/*_R1_001.fastq > samples_raw_fastqc.txt

# Submit FastQC job
sbatch 01_fastqc.sh
```

#### 1.2 Trim and Process Reads
```bash
# Create sample list (update array number in script based on sample count)
ls raw_files_combined/*_R1_001.fastq | sed 's/_R1_001.fastq//' | xargs -n1 basename > samples.txt

# Update SLURM --array parameter in 02_processing_reads_paired.sh
# Set: --array=0-N where N = (number of samples - 1)

sbatch 02_processing_reads_paired.sh
```
**Input:** `raw_files_combined/*_R[1,2]_001.fastq`  
**Output:** `trimmed_reads/*_R[1,2]_paired.fastq`, `trimmed_reads/*_R[1,2]_unpaired.fastq`

#### 1.3 FastQC on Trimmed Reads
```bash
# Create sample list for trimmed reads
ls trimmed_reads/*_R1_paired.fastq > samples_trimmed_fastqc.txt

sbatch 01_fastqc.sh
```

**Output:** `trimmed_fastqc/` directory with FastQC reports

#### 1.4 MultiQC Summary
```bash
sbatch 01b_multiqc.sh
```

---

### **PHASE 2: Contamination & Vector Removal**

#### 2.1 Remove Vector Contamination
```bash
# Update SLURM --array parameter based on sample count
sbatch 03_remove_vector_contamination_paired.sh
```
**Input:** `trimmed_reads/*_R[1,2]_paired.fastq`  
**Process:** 
  - Maps reads to UniVec_Core database using BWA
  - Filters out contaminating sequences using BLAT or custom script
  - Uses `1_BLAT_Filter_DB_paired.py` or `2_Infernal_Filter_paired.py`

**Output:** `quality_vector_rm_reads/*` (clean reads), `quality_vector_rm_reads/*_contam` (contaminated reads)

#### 2.2 Remove Host Reads (Optional)
download and index host genome first, change the location of the database in the script

```bash
sbatch 03b_remove_host_reads_paired.sh
```

---

### **PHASE 3: rRNA Removal & mRNA Extraction**

#### 3.1 Run rRNA Filtering with Infernal
```bash
# Update SLURM --array parameter (typically --array=1-N for 1-based indexing)
sbatch 04_rRNA_removal_array.sh
```
**Input:** `quality_vector_rm_reads/*`  
**Process:** 
  - Uses Infernal + OMM12 rRNA database
  - Searches for rRNA sequences in reads
  - Separates rRNA from mRNA

**Output:** `rRNA_rm_OMM12_reads/*_mRNA.fastq`, `rRNA_rm_OMM12_reads/*_rRNA.fastq`

#### 3.2 Final rRNA Filtering & Quality Control
```bash
sbatch 05_rRNA_filter_paired.sh
```
**Output:** `rRNA_filtered/*_R[1,2]_unique_mRNA.fastq`, `rRNA_filtered/*_R[1,2]_unique_rRNA.fastq`

#### 3.3 Visualize rRNA Removal Statistics
```bash
Rscript 05b_create_rRNA_overview_graphs.R
```
**Output:** Graphs showing read loss during rRNA removal

---

### **PHASE 4: Read Annotation & Quantification**

#### 4.1 Annotate mRNA Reads to Pseudogenome
```bash
# Create index of R1 files
ls rRNA_filtered/*_R1_unique_mRNA.fastq | wc -l  # Check count

# Update SLURM --array parameter
sbatch 06_annotate_mRNA_paired.sh
```
**Input:** `rRNA_filtered/*_R[1,2]_unique_mRNA.fastq`  
**Output:** `annotation_results/{SAMPLE}/*_annotation_bwa.sam`

#### 4.2 Annotate rRNA Reads (Optional)
```bash
sbatch 06b_annotate_rRNA_paired.sh
```

#### 4.3 Count Reads per Gene
```bash
sbatch 07_count_reads_paired.sh
```
**Process:** Uses featureCounts with GFF3 annotations  
**Output:** `annotation_results/{SAMPLE}/gene_counts.txt`

---

### **PHASE 5: Statistical Analysis & Functional Annotation **

#### 5.1 Merge Read Counts Across Samples
```bash
Rscript 08_merge_counts.R
```
**Input:** `annotation_results/{SAMPLE}/gene_counts.txt`  
**Output:** 
  - `results/gene_counts_wide_mRNA.csv` (sample columns, gene rows)
  - `results/gene_counts_wide_metadata.csv` (with sample metadata)

#### 5.2 Principal Component Analysis (PCA)
```bash
Rscript 09_PCA.R
```
**Input:** `results/gene_counts_wide_mRNA.csv`  
**Output:** 
  - PCA plots
  - `results/gene_counts_wide_metadata.csv` (with PCA scores)

#### 5.3 Differential Expression Analysis (DESeq2)
```bash
Rscript 10_DESeq2_taxon_normalised.R
```
**Input:** `results/gene_counts_wide_mRNA.csv`  
**Output:** 
  - DESeq2 results files
  - Normalized counts
  - MA plots, heatmaps
  
**Note:** This script implements taxon-normalized DESeq2 approach (Klingenberg & Meinicke 2017)

---

### **PHASE 6: Functional & Pathway Analysis**

#### 6.1 Annotate Genomes with Functional Information
```bash
Rscript 11a_functional_annotation_genomes.R
```

#### 6.1.1 Functional Annotation with CAYMAN
```bash
# Requires:
# - .faa protein files in "cleaned_faa files/" directory
# - CAYMAN database in "v3/" directory
# - Cutoffs file: "v3/cutoffs.csv"

bash run_cayman.sh
```

**Output:** `cleaned_faa files/cayman/*_cayman.csv` (functional annotations CAZymes)

#### 6.1.2 Functional Annotation with eggnog
```bash
# Requires:
# - .faa protein files in "cleaned_faa files/" directory
# - eggnog database in cayman directory

bash eggnog_map_loop.sh
```
**Output:** `eggnog_annotations/*.emapper.annotations` (functional annotations KEGG)

#### 6.2 Merge Functional Annotations with DESeq2 Results
```bash
Rscript 11b_merge_functional_annotations_DESEq2.R
```

#### 6.3 Extract CAZyme Protein Sequences
```bash
Rscript 11c_extract_cazyme_protein_sequences.R
```

#### 6.4 Volcano Plots
```bash
Rscript 12_volcanoplots.R
```
**Output:** Volcano plots highlighting significantly DEGs

#### 6.5 KEGG Pathway Enrichment
```bash
Rscript 13a_KEGG_pathway_enrichment.R
```

#### 6.6 KEGG Pathway Visualization
```bash
Rscript 13b_KEGG_pathview_maps.R
```
**Output:** Annotated KEGG pathway maps

#### 6.7 Gene Set Enrichment Analysis (GSEA)
```bash
Rscript 13c_GSEA_enrichment.R
```

---

## Output Structure

After completing the full pipeline, you'll have:

```
results/
├── gene_counts_wide_mRNA.csv          # Raw counts matrix
├── gene_counts_wide_metadata.csv      # Counts with sample metadata
├── pca_plots.pdf                      # PCA visualizations
├── DESeq2_results.csv                 # DE results
├── volcano_plots.pdf                  # Volcano plots
├── pathway_enrichment/                # KEGG enrichment results
├── pathview_maps/                     # KEGG pathway visualizations
└── gsea_results/                      # Gene set enrichment results

annotation_results/
└── {SAMPLE}/
    ├── gene_counts.txt                # featureCounts output
    ├── {SAMPLE}_annotation_bwa.sam    # BWA mapping results
    └── ...

trimmed_fastqc/                         # Quality reports
```

---

## Key Dependencies

- **Bioinformatics Tools:**
  - Trimmomatic (read trimming)
  - FastQC / MultiQC (quality control)
  - BWA (read mapping)
  - Samtools (BAM/SAM manipulation)
  - VSEARCH (homology search)
  - Infernal (rRNA detection)
  - featureCounts (read counting)
  - CAYMAN (functional annotation CAZymes)
  - eggnog (functional annotation KEGG)

- **R Packages:**
  - DESeq2
  - ggplot2
  - dplyr
  - tidyverse
  - pheatmap
  - pathview
  - fgsea (GSEA)

- **Python:**
  - BioPython
  - numpy, pandas

---

## Tips for New Users

1. **Start with Phase 0:** Always prepare databases, pseudogenome, and sample files first
2. **Array Job Sizing:** Update `--array` parameter based on your sample count
3. **Resource Allocation:** Adjust `--time`, `--mem-per-cpu`, and `--cpus-per-task` based on your HPC environment
4. **Check Logs:** Monitor `.log` files after each SLURM submission
5. **Sample Names:** Keep sample naming consistent throughout (e.g., `Sample_S01`, not `Sample-S01`)
6. **Test Run:** Run one sample through each step to validate before submitting full array jobs
7. **Storage:** Intermediate files can be large (~1-10 GB per sample); ensure adequate disk space
8. **Customize R Scripts:** Edit metadata extraction (media, replicate) in R scripts to match your experimental design

---

## References

- Klingenberg & Meinicke (2017): "How to normalize metatranscriptomic count data for differential expression analysis"
- Infernal: <http://eddylab.org/infernal/>
- DESeq2: <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>
- BWA: <http://bio-bwa.sourceforge.net/>
