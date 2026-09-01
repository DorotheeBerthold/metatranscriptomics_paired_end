# ==============================================================================
# Functional Annotation Integration Pipeline
# Integration of KAAS, KofamKOALA, eggNOG-mapper, dbCAN, and Cayman annotations
# ==============================================================================

library(tidyverse)
library(fs)

# 1. Directory & Parameter Setup -----------------------------------------------

org_list <- c("I46", "I48", "YL2", "YL31", "YL32", "YL44", "YL45", "YL58", "KB18", "KB1")

# Base directories
dir_gbff       <- file.path("OMM_files", "OMM_gbff")
dir_temp_anno  <- file.path("OMM_files", "temp_annotations")
dir_fasta      <- file.path("OMM_files", "fasta_protein")
dir_out        <- "annotated_genomes"

# Create output directories if they do not exist
dir_create(c(dir_temp_anno, dir_fasta, dir_eggnog, dir_out))

# 2. Extract CDS from GBFF & Generate Clean FASTAs -----------------------------

parse_cds_block <- function(block) {
  block_clean <- str_squish(block)
  
  locus_tag   <- str_match(block_clean, "/locus_tag=\"(.*?)\"")[, 2]
  gene        <- str_match(block_clean, "/gene=\"?([^\";\\s]+)")[, 2]
  protein_id  <- str_match(block_clean, "/protein_id=\"(.*?)\"")[, 2]
  
  translation <- str_match(block, "(?s)/translation=\"(.*?)\"")[, 2] %>%
    str_remove_all("\\s+")
  
  tibble(locus_tag, gene, protein_id, translation)
}

for (org in org_list) {
  message("Extracting proteins for: ", org)
  
  gbff_file <- file.path(dir_gbff, paste0(org, ".gb"))
  if (!file_exists(gbff_file)) next
  
  gbff_text <- read_file(gbff_file)
  feature_blocks <- str_split(gbff_text, "\n(?=     [A-Za-z]+)")[[1]]
  cds_blocks     <- feature_blocks[str_starts(feature_blocks, "     CDS")]
  
  anno <- map_dfr(cds_blocks, parse_cds_block) %>%
    filter(!is.na(protein_id), !is.na(translation))
  
  # Write raw and cleaned FASTA files directly in R
  fasta_lines <- paste0(">", anno$protein_id, "\n", anno$translation)
  write_lines(fasta_lines, file.path(dir_fasta, paste0(org, "_cleaned.faa")))
  
  # Save CDS reference table
  write_rds(select(anno, -translation), file.path(dir_temp_anno, paste0(org, "_anno.rds")))
}

# 3. Combine Multi-Tool Annotations --------------------------------------------

for (org in org_list) {
  message("Combining annotations for: ", org)
  
  anno_file <- file.path(dir_temp_anno, paste0(org, "_anno.rds"))
  if (!file_exists(anno_file)) next
  anno <- read_rds(anno_file)
  
  # --- KAAS ---
  kaas_map <- read.delim(paste0("OMM_files/fasta_protein/", org, "_kaas.txt"), 
                         header = FALSE, sep = "\t",
                         col.names = c("protein_id", "kegg_id_kegg"))
  
  # remove empty lines
  kaas_map <- kaas_map[kaas_map$kegg_id_kegg != "", ]
  
  # --- KofamKOALA ---
  koala_map <- read.delim(paste0("OMM_files/fasta_protein/", org, "_koala.txt"), header = FALSE, sep = "\t",
                              col.names = c("protein_id", "kegg_id_koala"))
  
  # remove empty lines
  koala_map <- koala_map[koala_map$kegg_id_koala != "", ]
  
  # --- eggNOG-mapper ---
  lines <- readLines(paste0("OMM_files/eggnog_annotations2/", org, ".emapper.annotations"))
  header_line_index <- which(grepl("^#query", lines))[1]
  
  eggnog <- read.delim(paste0("OMM_files/eggnog_annotations2/", org, ".emapper.annotations"),
                       header = TRUE, sep = "\t", quote = "", fill = TRUE,
                       skip = header_line_index - 1, stringsAsFactors = FALSE)
  
  eggnog_map <- eggnog[, c("X.query", "KEGG_ko")]
  colnames(eggnog_map) <- c("protein_id","kegg_id_eggnog")
  eggnog_map$kegg_id_eggnog <- gsub("ko:", "", eggnog_map$kegg_id_eggnog)
  
  eggnog_map <- eggnog_map[!grepl("^#", eggnog_map$protein_id), ]
  
  eggnog_map <- eggnog_map %>%
    filter(kegg_id_eggnog != "-" & !is.na(kegg_id_eggnog))
  
  # --- dbCAN ---
  dbcan <- read.delim(paste0("OMM_files/fasta_protein/", org, "_dbcan.txt"), header = TRUE, sep = "\t")
  
  dbcan_map <- dbcan[, c("Gene.ID", "HMMER")]
  colnames(dbcan_map) <- c("protein_id", "dbcan_family")
  
  # keep only family
  dbcan_map$dbcan_family <- sub("\\(.*\\)$", "", dbcan_map$dbcan_family)
  
  # --- Cayman ---
  cayman_map <- read.csv(paste0("OMM_files/fasta_protein/", org, "_cleaned_cayman.csv"))
  cayman_map <- cayman_map[, c("sequenceID", "family")]
  colnames(cayman_map) <- c("protein_id", "cayman_family")
  
  # --- Merge & Consensus Logic ---
  merged_ko <- anno %>%
    select(protein_id) %>%
    left_join(kaas_map, by = "protein_id") %>%
    left_join(eggnog_map, by = "protein_id") %>%
    left_join(koala_map, by = "protein_id") %>%
    left_join(dbcan_map, by = "protein_id") %>%
    left_join(cayman_map, by = "protein_id") %>%
    mutate(
      kegg_consensus = case_when(
        # 1. Majority consensus (2 out of 3 agree)
        !is.na(kegg_id_koala) & kegg_id_koala == kegg_id_eggnog ~ kegg_id_koala,
        !is.na(kegg_id_koala) & kegg_id_koala == kegg_id_kegg   ~ kegg_id_koala,
        !is.na(kegg_id_eggnog) & kegg_id_eggnog == kegg_id_kegg ~ kegg_id_eggnog,
        
        # 2. Priority fallback (KofamKOALA > eggNOG v3 > KAAS)
        !is.na(kegg_id_koala)  ~ kegg_id_koala,
        !is.na(kegg_id_eggnog) ~ kegg_id_eggnog,
        !is.na(kegg_id_kegg)   ~ kegg_id_kegg,
        
        TRUE ~ NA_character_
      )
    )
  
  # Output complete annotation table
  final_annot <- left_join(anno, merged_ko, by = "protein_id")
  write_csv(final_annot, file.path(dir_out, paste0(org, "_annotated.csv")))
}
