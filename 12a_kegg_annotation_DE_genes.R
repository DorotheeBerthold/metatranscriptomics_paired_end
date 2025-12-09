# KEGG extraction =======================================================================
# Using KAAS and eggNOG to annotate DE genes per organism with KEGG Orthology (KO) IDs
# ========================================================================================


# Libraries ------------------------------------------------------------------------------

library(tidyverse)
library(stringr)

# Data import ----------------------------------------------------------------------------

# Import shrunken results from DESEq2 pipeline
results <- read.csv("results/DB092_deseq2_results_shrunken.csv", row.names = 1)
results$genes <- rownames(results)

# Filter foreach organism by using grepl on rownames
org_list <- c("I46","I48","YL2","YL31", "YL32", "YL44", "YL45", "YL58")

# separate genes column into two based on .
bacteria_results <- results %>% 
  separate(genes, into = c("gene", "bacteria"), sep = "\\.") %>% 
  filter(!is.na(padj)) %>%
  mutate(diffexpressed = recode(diffexpressed, 
                                "NO" = "", 
                                "UP" = "+fiber", 
                                "DOWN" = "-fiber")) 

sig_results <- bacteria_results %>% 
  filter(diffexpressed != "")

# Generate FASTA files for KAAS & eggNOG submission -----------------------------------

org_list <- c("I46","I48","YL2","YL31", "YL32", "YL44", "YL45", "YL58")

for(org in org_list) {
  message("Preparing FASTA for organism: ", org)
  
  # Filter results for this organism
  org_results <- sig_results %>% filter(bacteria == org)
  
  # Read GBFF file
  gbff_lines <- readLines(paste0("OMM_gbff/", org, ".gbff"))
  gbff_text <- paste(gbff_lines, collapse = "\n")
  
  # Split into CDS blocks
  cds_blocks <- str_split(gbff_text, "     CDS")[[1]][-1]
  
  # Parse CDS
  parse_cds <- function(block){
    block <- str_squish(block)
    locus_tag  <- str_match(block, "/locus_tag=\"(.*?)\"")[,2]
    protein_id <- str_match(block, "/protein_id=\"(.*?)\"")[,2]
    translation <- str_match(block, "(?s)/translation=\"(.*?)\"")[,2] %>%
      str_replace_all("\n", "") %>% str_squish()
    tibble(locus_tag, protein_id, translation)
  }
  
  anno <- map_df(cds_blocks, parse_cds) %>%
    filter(!is.na(protein_id), !is.na(translation))
  
  # Write FASTA
  fasta_lines <- paste0(">", anno$protein_id, "\n", anno$translation)
  writeLines(fasta_lines, paste0("fasta_protein/", org, "_proteins.faa"))
  
  message("FASTA written for ", org, ". Submit to KAAS and eggNOG before continuing.")
  
  # Save the annotation without translation for later
  write_rds(anno %>% select(-translation), paste0("temp_annotations/", org, "_anno.rds"))
  
  # Break here to allow manual submission
  # Uncomment the next line to pause after each organism
  # readline(prompt = "Press [Enter] after submitting to KAAS/eggnog to continue...")
}

# Submit to KAAS & eggNOG ----------------------------------------------------------------
# 1. KAAS: https://www.genome.jp/tools/kaas/
#    - Upload the generated FASTA files
#    - Select "BBH" method and appropriate taxonomic group
#    - Download the results as "<organism>_kaas.txt" and save in "fasta_protein/" folder

# 2. eggNOG-mapper: https://emapper.embl.de/
#    - Use 12b_eggnog_map_loop.sh script to batch submit FASTA files
#    - Download the results as "<organism>.emapper.annotations" and save in "eggnog_annotations/" folder

# Combine annotations with DE results  -------------------------------------------------


for(org in org_list) {
  message("Processing annotations for organism: ", org)
  
  # Read original results & saved annotations
  org_results <- sig_results %>% filter(bacteria == org)
  anno <- read_rds(paste0("temp_annotations/", org, "_anno.rds"))
  
  org_results_annot <- inner_join(org_results, anno, by = c("gene" = "locus_tag"))
  
  # Read KAAS mapping
  kegg_map <- read.delim(paste0("fasta_protein/", org, "_kaas.txt"), 
                         header = FALSE, sep = "\t",
                         col.names = c("protein_id", "kegg_id_kegg"))
  
  # Read eggNOG mapping
  lines <- readLines(paste0("eggnog_annotations/", org, ".emapper.annotations"))
  header_line_index <- which(grepl("^#query", lines))[1]
  
  eggnog <- read.delim(paste0("eggnog_annotations/", org, ".emapper.annotations"),
                       header = TRUE, sep = "\t", quote = "", fill = TRUE,
                       skip = header_line_index - 1, stringsAsFactors = FALSE)
  
  eggnog_map <- eggnog[, c("X.query", "KEGG_ko")]
  colnames(eggnog_map) <- c("protein_id","kegg_id_eggnog")
  eggnog_map$kegg_id_eggnog <- gsub("ko:", "", eggnog_map$kegg_id_eggnog)
  
  # Combine KO mappings
  ko_map <- full_join(kegg_map, eggnog_map, by = "protein_id") %>%
    mutate(kegg_id = ifelse(kegg_id_kegg != "", kegg_id_kegg, kegg_id_eggnog)) %>%
    select(protein_id, kegg_id)
  
  org_results_annot_ko <- left_join(org_results_annot, ko_map, by = "protein_id")
  
  # Filter proteins with KO
  org_results_final <- org_results_annot_ko %>% filter(grepl("K", kegg_id))
  
  # Write output
  write.csv(org_results_final, paste0("results/DE_genes_", org, "_annotated_KOs.csv"), row.names = FALSE)
}
