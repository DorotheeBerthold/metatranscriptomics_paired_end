# ==============================================================================
# Automated KEGG Pathview Visualization Pipeline
# Generates pathway map diagrams for up and down regulated datasets
# ==============================================================================

library(pathview)
library(stringr)
library(dplyr)

# 1. Helper Function to Generate Pathview Maps --------------------------------

generate_pathview_maps <- function(df, output_dir, fill_color = "red") {
  if (is.null(df) || nrow(df) == 0) {
    message("No entries to process for: ", output_dir)
    return(NULL)
  }
  
  # Ensure target directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Guarantee working directory restoration even if an error occurs
  original_wd <- getwd()
  on.exit(setwd(original_wd), add = TRUE)
  
  setwd(output_dir)
  message("Generating Pathview maps in: ", getwd())
  
  for (i in 1:nrow(df)) {
    current_row <- df[i, ]
    
    # Clean IDs and strain labels
    numeric_pathway_id <- str_remove(current_row$ID, "ko")
    safe_bacteria      <- str_extract(current_row$Bacteria, "^[A-Za-z0-9]+")
    
    # Parse K-numbers into named numeric vector
    gene_list <- unlist(strsplit(current_row$geneID, "/"))
    plot_data <- setNames(rep(1, length(gene_list)), gene_list)
    
    message("  -> Pathway: ", current_row$ID, " | Strain: ", safe_bacteria)
    
    # Render KEGG native map
    tryCatch({
      pathview(
        gene.data   = plot_data,
        pathway.id  = numeric_pathway_id,
        species     = "ko",
        kegg.native = TRUE,
        out.suffix  = safe_bacteria,
        limit       = list(gene = 1),
        low         = "white",
        high        = fill_color
      )
    }, error = function(e) {
      message("     [Error] Failed to generate map: ", e$message)
    })
  }
}

# 2. Automated Batch Processing ------------------------------------------------

# Process +Fiber (Upregulated)
generate_pathview_maps(
  df         = master_up_df,
  output_dir = "plots/pathview_maps_up2",
  fill_color = "red"
)

# Process -Fiber (Downregulated)
generate_pathview_maps(
  df         = master_down_df,
  output_dir = "plots/pathview_maps_down2",
  fill_color = "blue"
)
