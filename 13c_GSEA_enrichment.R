# ==============================================================================
# Functional Pathway Enrichment Analysis (GSEA) Across Bacterial Strains
# Automated GSEA execution, ridgeplot generation, and master result compilation
# ==============================================================================

library(clusterProfiler)
library(enrichplot)
library(tidyverse)

# 1. Setup & Directory Initialization ------------------------------------------

dir_results <- "results"
dir_plots   <- "plots"

dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

file_de_results <- file.path(dir_results, "DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed2.csv")
bacteria_results <- read_csv(file_de_results, show_col_types = FALSE)

org_list <- c("I46", "I48", "YL2", "YL31", "YL32", "YL44", "YL45", "YL58")

# 2. Main GSEA Execution Loop --------------------------------------------------

gsea_accumulator <- list()

for (org in org_list) {
  message("==================================================")
  message("Processing GSEA for organism: ", org)
  
  # Filter dataset for target organism
  org_data <- bacteria_results %>%
    filter(str_detect(bacteria, org))
  
  if (nrow(org_data) == 0) {
    message("  -> WARNING: No rows found matching '", org, "'. Skipping.")
    next
  }
  
  # Filter valid KEGG annotations
  org_data_clean <- org_data %>%
    filter(!is.na(kegg_consensus), kegg_consensus != "", kegg_consensus != "<NA>")
  
  if (nrow(org_data_clean) < 15) {
    message("  -> Skipping: Insufficient mapped KEGG terms (", nrow(org_data_clean), " genes) for GSEA.")
    next
  }
  
  # Prepare ranked named numeric vector (Log2FC in descending order)
  gsea_input <- org_data_clean %>%
    group_by(kegg_consensus) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(log2FoldChange)) %>%
    select(kegg_consensus, log2FoldChange) %>%
    deframe()
  
  # Run GSEA execution wrapped in tryCatch
  gsea_res <- tryCatch({
    gseKEGG(
      geneList     = gsea_input,
      organism     = "ko",
      keyType      = "kegg",
      pvalueCutoff = 0.05,
      minGSSize    = 10,
      maxGSSize    = 500,
      verbose      = FALSE
    )
  }, error = function(e) {
    message("  -> GSEA error for ", org, ": ", e$message)
    return(NULL)
  })
  
  # Process and save results if significant pathways found
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    
    # Prepend ID to Description to prevent factor level collision during plotting
    gsea_res@result$Description <- paste0(gsea_res@result$ID, ": ", gsea_res@result$Description)
    
    full_name <- unique(org_data$bacteria)[1]
    
    gsea_df <- as.data.frame(gsea_res) %>%
      mutate(
        org_code = org,
        bacteria = full_name
      ) %>%
      select(org_code, bacteria, everything())
    
    gsea_accumulator[[org]] <- gsea_df
    
    # Save individual strain CSV
    write_csv(gsea_df, file.path(dir_results, paste0("GSEA_results_", org, "2.csv")))
    
    # Render and save Ridgeplot
    p_ridge <- ridgeplot(gsea_res, showCategory = 12) +
      labs(
        x = expression(Log[2]~"Fold Change Rank Distribution"),
        title = paste("Pathway Expression Shifts:", org),
        subtitle = full_name
      ) +
      theme_minimal() +
      theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(face = "italic", size = 11, color = "grey30")
      )
    
    ggsave(
      filename = file.path(dir_plots, paste0("GSEA_ridgeplot_", org, "_2.png")),
      plot     = p_ridge,
      width    = 8,
      height   = 6
    )
    
    message("  -> SUCCESS: Saved GSEA table and ridgeplot for ", org)
    
  } else {
    message("  -> No significantly enriched pathways (p.adjust < 0.05) found for ", org)
  }
}

# 3. Master Dataset Consolidation ----------------------------------------------

master_gsea_df <- bind_rows(gsea_accumulator)
write_csv(master_gsea_df, file.path(dir_results, "GSEA_master_results_all_organisms2.csv"))
