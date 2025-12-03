# Merge gene_counts from all different samples

{
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
}

# List all gene_counts.txt files recursively
files <- list.files("annotation_results", pattern = "gene_counts.txt", full.names = TRUE, recursive = TRUE)

# Function to read a file and add sample_id based on folder name
read_with_sample <- function(file_path) {
  sample_id <- basename(dirname(file_path))
  
  df <- read_delim(file_path, skip = 1, show_col_types = F) 
  colnames(df)[7] <- sample_id
  
  df <- df %>%
    select(Geneid, count = all_of(sample_id), Chr) %>%
    mutate(sample_id = sample_id)
  return(df)
}

file1 <- read_delim("annotation_results/APFmod-MC-1_S01/gene_counts.txt", skip = 1, show_col_types = F)
#colnames(file1)[7] <- "APFmod-MC-1_S01"

# Filter out all Geneid with start 1
full_genome <- file1 %>%
  filter(Start == "1") 
# doesn't work as I48 genome does not start at 1

# Read all files and combine
merged_long <- as.data.frame(map_dfr(files, read_with_sample))
merged_long$gene_id <- paste0(merged_long$Geneid, "-", merged_long$Chr)

# Convert to wide format
merged_wide <- merged_long %>%
  pivot_wider(names_from = sample_id, values_from = count)

merged_wide <- as.data.frame(merged_wide)
rownames(merged_wide) <- merged_wide$gene_id

merged_wide <- merged_wide %>% 
  select(-c(Geneid, gene_id, Chr))

write.csv(merged_wide, file= "results/gene_counts_wide_mRNA.csv")

# remove rows where all columns have 0
merged_wide_filtered <- merged_wide[rowSums(merged_wide != 0) > 0, ]

write.csv(merged_wide_filtered, "results/gene_counts_wide_filtered.csv")
                                    