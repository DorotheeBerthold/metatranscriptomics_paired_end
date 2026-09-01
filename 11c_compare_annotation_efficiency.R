library(tidyverse)

# 1. Load Datasets -------------------------------------------------------------
f1_path <- "results/DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed.csv"
f2_path <- "results/DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed2.csv"

df1 <- read_csv(f1_path, show_col_types = FALSE)
df2 <- read_csv(f2_path, show_col_types = FALSE)

# 2. Summarize Annotation Function ---------------------------------------------
summarize_tools <- function(df, label) {
  # Columns to check for annotations
  target_cols <- intersect(
    c("kegg_id_kegg", "kegg_id_eggnog", "kegg_id_koala", "kegg_consensus", "dbcan_family", "cayman_family"),
    colnames(df)
  )
  
  df %>%
    group_by(bacteria) %>%
    summarize(
      total_genes = n(),
      across(
        all_of(target_cols),
        ~ sum(!is.na(.) & . != "" & . != "-" & . != "NA"),
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = -c(bacteria, total_genes),
      names_to = "tool",
      values_to = paste0("count_", label)
    )
}

# 3. Calculate Counts & Merge Comparison ---------------------------------------
stats1 <- summarize_tools(df1, "file1")
stats2 <- summarize_tools(df2, "file2")

comparison <- full_join(stats1, stats2, by = c("bacteria", "tool")) %>%
  mutate(
    count_file1 = replace_na(count_file1, 0),
    count_file2 = replace_na(count_file2, 0),
    total_genes = coalesce(total_genes.y, total_genes.x),
    pct_file1   = round((count_file1 / total_genes) * 100, 1),
    pct_file2   = round((count_file2 / total_genes) * 100, 1),
    diff_count  = count_file2 - count_file1
  ) %>%
  select(bacteria, tool, total_genes, count_file1, pct_file1, count_file2, pct_file2, diff_count)

# 4. Reshape to Wide Matrix (Organism x Tool) ----------------------------------
wide_summary <- comparison %>%
  pivot_wider(
    id_cols = c(bacteria, total_genes),
    names_from = tool,
    values_from = c(count_file2, pct_file2)
  )

# 5. Print Results & Export ----------------------------------------------------
print(comparison, n = Inf)

write_csv(comparison, "annotation_comparison_long.csv")
write_csv(wide_summary, "annotation_comparison_wide.csv")