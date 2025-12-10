# KEGG enrichment  =======================================================================
# Generate KEGG pathway enrichment plots for multiple organisms
# ========================================================================================

# Libraries --------------------------------------------------------------------------------

library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(purrr)

# Functions ------------------------------------------------------------------------------
split_kegg <- function(x) {
  ids <- unlist(strsplit(as.character(x), ","))
  ids <- trimws(ids)            # remove leading/trailing spaces
  ids <- ids[ids != ""]         # drop empty entries
  ids <- ids[!is.na(ids)]       # drop NAs
  unique(ids)
}

# Data import & loop ---------------------------------------------------------------------
org_list <- c("I46","I48","YL2","YL31", "YL32", "YL44", "YL45", "YL58")

for (org in org_list) {
  
  message("Processing: ", org)
  
  # -------------------------------
  # Load file
  # -------------------------------
  file_path <- paste0("results/DE_genes_", org, "_annotated_KOs.csv")
  org_results_final <- read.csv(file_path)
  
  # -------------------------------
  # Split by regulation
  # -------------------------------
  upregulated <- org_results_final %>% 
    filter(diffexpressed == "+fiber")
  
  downregulated <- org_results_final %>% 
    filter(diffexpressed == "-fiber")
  
  # Background universe (optional)
  # background_genes <- org_results_final %>% filter(diffexpressed == "")
  # background_kegg <- split_kegg(background_genes$kegg_id)
  
  # -------------------------------
  # Extract KEGG IDs
  # -------------------------------
  up_kegg <- split_kegg(upregulated$kegg_id)
  down_kegg <- split_kegg(downregulated$kegg_id)
  
  # -------------------------------
  # KEGG enrichment
  # -------------------------------
  enrich_up <- enrichKEGG(
    gene = up_kegg,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = 0.05
    # universe = background_kegg
  )
  
  enrich_down <- enrichKEGG(
    gene = down_kegg,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = 0.05
    # universe = background_kegg
  )
  
  # -------------------------------
  # Convert to df
  # -------------------------------
  up_df <- enrich_up@result %>% 
    mutate(direction = "+ Fiber")
  
  down_df <- enrich_down@result %>% 
    mutate(direction = "- Fiber")
  
  # -------------------------------
  # Combine
  # -------------------------------
  enrich_combined <- bind_rows(up_df, down_df)
  
  # Keep top 20 per group
  enrich_combined <- enrich_combined %>%
    group_by(direction) %>%
    slice_min(pvalue, n = 20) %>%
    ungroup()
  
  # -------------------------------
  # Plot
  # -------------------------------
  p <- ggplot(enrich_combined,
              aes(x = Count,
                  y = reorder(Description, -pvalue),
                  size = Count,
                  color = p.adjust)) +
    geom_point() +
    scale_color_continuous(trans = "log10", name = "Adjusted p-value") +
    facet_wrap(~ direction, scales = "free_y") +
    theme_bw() +
    labs(
      title = paste("KEGG Pathway Enrichment:", org),
      x = "Gene count",
      y = "Pathway"
    )
  
  # Save figure
  out_png <- paste0("plots/KEGG_pathway_enrichment_", org, ".png")
  ggsave(out_png, p, width = 10, height = 6)
  
  message("Finished: ", org, " → saved ", out_png)
}

# Create community level summary (optional) ---------------------------------------------

community_kegg_list <- lapply(org_list, function(org) {
  read.csv(paste0("results/DE_genes_", org, "_annotated_KOs.csv"))
})

# Combine all into one big dataframe
community_df <- dplyr::bind_rows(community_kegg_list)

upregulated <- community_df %>% 
  filter(diffexpressed == "+fiber")

downregulated <- community_df %>% 
  filter(diffexpressed == "-fiber")

up_kegg <- split_kegg(upregulated$kegg_id)
down_kegg <- split_kegg(downregulated$kegg_id)

# -------------------------------
# KEGG enrichment
# -------------------------------
enrich_up <- enrichKEGG(
  gene = up_kegg,
  organism = "ko",
  keyType = "kegg",
  pvalueCutoff = 0.05
  # universe = background_kegg
)

enrich_down <- enrichKEGG(
  gene = down_kegg,
  organism = "ko",
  keyType = "kegg",
  pvalueCutoff = 0.05
  # universe = background_kegg
)

# -------------------------------
# Convert to df
# -------------------------------
up_df <- enrich_up@result %>% 
  mutate(direction = "+ Fiber")

down_df <- enrich_down@result %>% 
  mutate(direction = "- Fiber")

# -------------------------------
# Combine
# -------------------------------
enrich_combined <- bind_rows(up_df, down_df)

# Keep top 20 per group
enrich_combined <- enrich_combined %>%
  group_by(direction) %>%
  slice_min(pvalue, n = 20) %>%
  ungroup()

ggplot(enrich_combined,
       aes(x = Count,
           y = reorder(Description, -pvalue),
           size = Count,
           color = p.adjust)) +
  geom_point() +
  scale_color_continuous(trans = "log10", name = "Adjusted p-value") +
  facet_wrap(~ direction, scales = "free_y") +
  theme_bw() +
  labs(
    title = paste("KEGG Pathway Enrichment Community"),
    x = "Gene count",
    y = "Pathway"
  )
# Save figure
ggsave("plots/KEGG_pathway_enrichment_community.png", width = 10, height = 6)

