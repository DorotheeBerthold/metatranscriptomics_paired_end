library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(ggfortify)


# For PCA: samples as colnames, metabolites as rownames
merged_wide_filtered <- read.csv("results/gene_counts_wide_mRNA.csv", row.names = 1)


# make all numeric
gene_matrix <- as.matrix(sapply(merged_wide_filtered, as.numeric))
rownames(gene_matrix) <- rownames(merged_wide_filtered)


# Replace NA values with 0
ms_results_data <- as.data.frame(t(gene_matrix))

ms_results_data$annotation <- rownames(ms_results_data)

ms_results_data <- ms_results_data %>%
  mutate(replicate = str_extract(annotation, "(?<=MC\\.)\\d+")) %>% 
  mutate(media = case_when(grepl("mod", annotation) ~"-fiber",
                           grepl("xylan", annotation) ~ "+fiber")) %>% 
  relocate(media, annotation, replicate) %>% 
    select(-annotation)
write.csv(ms_results_data, "results/gene_counts_wide_metadata.csv")

# make long
ms_results_long <- ms_results_data %>%
  pivot_longer(cols = -c(media, replicate),
               names_to = "gene",
               values_to = "count")

#write.csv(ms_results_long, "results/gene_counts_long.csv")

pca_data <- t(gene_matrix)
pca_data[is.na(pca_data)] <- 0
# Remove columns with zero variance
pca_data <- pca_data[, apply(pca_data, 2, var) != 0]
pca_result <- prcomp(pca_data,center = TRUE,  scale. = TRUE, )
screeplot(pca_result)

# Extracting scores
pca_scores <- as.data.frame(pca_result$x)

# Add metadata
pca_scores$media <- as.factor(ms_results_data$media)
pca_scores$replicate <- as.factor(ms_results_data$replicate)


autoplot(pca_result, data = pca_scores, colour = "media", size = 7, loadings = F) + 
  theme_classic() +
  #scale_shape_manual(values = custom_shapes) +
  labs(title = "Genes in the different diets")


ggsave("plots/pca_genes_different_diets.png", width = 6, height = 5)
