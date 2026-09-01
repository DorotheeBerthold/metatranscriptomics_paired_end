

# Import shrunken results from DESEq2 pipeline
results <- read.csv("results/DB092_deseq2_results_shrunken.csv", row.names = 1)
# Annotate scaled matrix
#results <- read.csv("results/DESeq2_taxon_specific_scaled_matrix.csv", row.names = 1)
results$genes <- rownames(results)

# separate genes column into two based on .
results <- results %>% 
  separate(genes, into = c("locus_tag", "bacteria"), sep = "\\.")


# Read in annotated genomes
files <- list.files("annotated_genomes/", pattern = "*_annotated.csv", full.names = TRUE)

# read them all in and rbind
annotated_bacteria <- map_df(files, read.csv)

write.csv(annotated_bacteria, "results/annotated_bacteria_all.csv", row.names = F)


bacteria_results_annotated <- inner_join(results, annotated_bacteria, by = c("locus_tag"))

bacteria_results_annotated <- bacteria_results_annotated %>%
  mutate(bacteria = recode(bacteria,
                           
                           "I48" = "Bacteroides caecimuris I48",
                           "YL58" = "Blautia coccoides YL58",
                           "YL32" = "Enterocloster clostridioformis YL32",
                           "YL45" = "Turicimonas muris YL45",
                           "I46" = "Clostridium innocuum I46",
                           "YL44" = "Akkermansia muciniphila YL44",
                           "YL31" = "Flavonifractor plautii YL31",
                           "YL2" = "Bifidobacterium animalis YL2"
  )) %>% 
  mutate(diffexpressed = if_else(diffexpressed == "UP", "+fiber", 
                                 if_else(diffexpressed == "DOWN", "-fiber", diffexpressed)))


write.csv(bacteria_results_annotated, "results/DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed2.csv", row.names = FALSE)
