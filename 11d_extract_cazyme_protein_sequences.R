# ==============================================================================
# CAZyme Domain Sequence Slicing Pipeline
# Extracts target CAZyme domain amino acid sequences based on Cayman coordinates
# ==============================================================================

library(Biostrings)
library(tidyverse)
library(fs)

# 1. Setup & Parameters -------------------------------------------------------

org_list <- c("I46", "I48", "YL2", "YL31", "YL32", "YL44", "YL45", "YL58", "KB18", "KB1")

dir_fasta <- file.path("OMM_files", "fasta_protein")
dir_out   <- file.path("OMM_files", "cazyme_slices")

# Automatically create output directory if missing
dir_create(dir_out)

# 2. Extraction Function -------------------------------------------------------

extract_cazy_slices <- function(org) {
  message("Extracting CAZyme slices for: ", org)
  
  cazy_file <- file.path(dir_fasta, paste0(org, "_cleaned_cayman.csv"))
  faa_file  <- file.path(dir_fasta, paste0(org, "_cleaned.faa"))
  
  # Safety check for missing files
  if (!file_exists(cazy_file) || !file_exists(faa_file)) {
    warning("Missing input files for ", org, ". Skipping...")
    return(NULL)
  }
  
  # Read Cayman domain output and FASTA file
  cazy <- read_csv(cazy_file, show_col_types = FALSE)
  faa  <- readAAStringSet(faa_file)
  
  # Create fast lookup reference frame
  seq_df <- tibble(
    sequenceID       = str_remove(names(faa), "\\s+.*$"), # Clean header titles
    protein_sequence = as.character(faa)
  )
  
  # Join and slice target domain sequences
  result <- cazy %>%
    mutate(sequenceID = str_remove_all(sequenceID, "^b'|'$")) %>%
    left_join(seq_df, by = "sequenceID") %>%
    mutate(
      # Vectorized slicing using str_sub eliminates slow rowwise execution
      slice    = str_sub(protein_sequence, start_protein, end_protein),
      organism = org
    )
  
  # Write individual organism CSV
  write_csv(result, file.path(dir_out, paste0(org, "_cazyme_slices.csv")))
  
  return(result)
}

# 3. Batch Execution & Consolidation -------------------------------------------

# Process all organisms and bind non-empty results
cazy_combined <- map(org_list, extract_cazy_slices) %>%
  compact() %>%
  bind_rows()

# Export combined dataset
write_csv(cazy_combined, file.path(dir_out, "all_cazyme_slices.csv"))

message("Finished extracting CAZyme domain slices across all organisms.")