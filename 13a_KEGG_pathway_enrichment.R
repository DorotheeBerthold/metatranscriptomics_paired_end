# ==============================================================================
# KEGG Pathway Enrichment Analysis & Cross-Organism Comparative Visualization
# Per-organism cluster comparison and community-wide pathway summary
# ==============================================================================

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(clusterProfiler)

# 1. Setup & Helper Functions --------------------------------------------------

# Ensure output directories exist
dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Parse comma-separated KEGG IDs into clean unique vectors
split_kegg <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x[!is.na(x)] %>%
    strsplit(",") %>%
    unlist() %>%
    trimws() %>%
    .[. != ""] %>%
    unique()
}

# Helper to process enrichResult objects into standardized data frames
parse_enrichment <- function(enrich_obj, org, gene_count) {
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) return(NULL)
  
  as.data.frame(enrich_obj) %>%
    filter(subcategory != "Global and overview maps") %>%
    mutate(
      Bacteria = paste0(org, "\n(n=", gene_count, ")"),
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
    )
}

# 2. Data Import & Configuration -----------------------------------------------

org_list <- c("I46", "I48", "YL2", "YL31", "YL32", "YL44", "YL45", "YL58")
file_de_results <- file.path("results", "DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed2.csv")

org_results <- read_csv(file_de_results, show_col_types = FALSE)

# Lists to accumulate community-wide enrichment tables
master_up_list   <- list()
master_down_list <- list()

# 3. Main Processing Loop ------------------------------------------------------

for (org in org_list) {
  
  message("Processing organism: ", org)
  
  # Filter organism-specific DE dataset
  org_results_final <- org_results %>% filter(str_detect(bacteria, org))
  
  # Extract gene subsets and KEGG KO lists
  upregulated   <- org_results_final %>% filter(diffexpressed == "+fiber")
  downregulated <- org_results_final %>% filter(diffexpressed == "-fiber")
  
  up_kegg   <- split_kegg(upregulated$kegg_consensus)
  down_kegg <- split_kegg(downregulated$kegg_consensus)
  bg_kegg   <- split_kegg(org_results_final$kegg_consensus)
  
  # --- Part A: Per-Organism CompareCluster Dotplots ---
  gene_list <- setNames(
    list(up_kegg, down_kegg),
    c(paste0("+ Fiber (Total: ", length(up_kegg), ")"), 
      paste0("- Fiber (Total: ", length(down_kegg), ")"))
  )
  
  kegg_comparison <- tryCatch({
    compareCluster(
      geneCluster  = gene_list,
      fun          = "enrichKEGG",
      organism     = "ko",
      keyType      = "kegg",
      universe     = bg_kegg,
      pvalueCutoff = 0.05
    )
  }, error = function(e) {
    message("  -> API/Enrichment error for ", org, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(kegg_comparison) && nrow(as.data.frame(kegg_comparison)) > 0) {
    # Remove broad global maps
    kegg_comparison@compareClusterResult <- kegg_comparison@compareClusterResult %>%
      filter(subcategory != "Global and overview maps")
    
    if (nrow(kegg_comparison@compareClusterResult) > 0) {
      p_ind <- dotplot(kegg_comparison, showCategory = 10, title = paste0("KEGG Pathways: ", org))
      ggsave(paste0("plots/KEGG_pathway_enrichment2_", org, ".png"), p_ind, width = 10, height = 6)
    }
  }
  
  # --- Part B: Accumulate Data for Community-Wide Plots ---
  
  # Upregulated (+Fiber)
  if (length(up_kegg) > 0) {
    enrich_up <- tryCatch({
      enrichKEGG(gene = up_kegg, organism = "ko", keyType = "kegg", universe = bg_kegg, pvalueCutoff = 0.05)
    }, error = function(e) NULL)
    
    master_up_list[[org]] <- parse_enrichment(enrich_up, org, length(up_kegg))
  }
  
  # Downregulated (-Fiber)
  if (length(down_kegg) > 0) {
    enrich_down <- tryCatch({
      enrichKEGG(gene = down_kegg, organism = "ko", keyType = "kegg", universe = bg_kegg, pvalueCutoff = 0.05)
    }, error = function(e) NULL)
    
    master_down_list[[org]] <- parse_enrichment(enrich_down, org, length(down_kegg))
  }
}

# 4. Master Datasets & Community Visualizations --------------------------------

master_up_df   <- bind_rows(master_up_list)
master_down_df <- bind_rows(master_down_list)

# Custom dotplot generator
plot_custom_dotplot <- function(df, condition_title) {
  if (is.null(df) || nrow(df) == 0) {
    message("No significant pathways found for ", condition_title)
    return(NULL)
  }
  
  ggplot(df, aes(x = Bacteria, y = Description)) +
    geom_point(aes(size = Count, fill = p.adjust), shape = 21, color = "black", stroke = 0.5) +
    scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +
    scale_size_continuous(range = c(3, 8)) +
    theme_bw() +
    labs(title = paste("KEGG Pathways Across Community:", condition_title), x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey92")
    )
}

# Plot community dotplots
p_up   <- plot_custom_dotplot(master_up_df, "+ Fiber")
p_down <- plot_custom_dotplot(master_down_df, "- Fiber")

if (!is.null(p_up))   ggsave("plots/KEGG_All_Bacteria_Up2.png", p_up, width = 10, height = 8)
if (!is.null(p_down)) ggsave("plots/KEGG_All_Bacteria_Down2.png", p_down, width = 10, height = 8)

# Save consolidated results
write_csv(master_up_df, "results/kegg_enrichment_all_bacteria_up2.csv")
write_csv(master_down_df, "results/kegg_enrichment_all_bacteria_down2.csv")
