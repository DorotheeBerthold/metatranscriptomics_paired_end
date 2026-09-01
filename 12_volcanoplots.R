# ==============================================================================
# Functional Categorization & Volcano Plot Visualization
# Integration of DE Results, Cayman Substrates, and KO Pathways
# ==============================================================================

library(tidyverse)
library(readxl)
library(ggrepel)
library(ggpubr)

# 1. File Paths & Parameter Setup ----------------------------------------------

file_de_results <- file.path("results", "DE_genes_all_bacteria_annotated_KOs_all_with_nodiffexpressed2.csv")
file_cayman     <- file.path("cayman", "cayman_annotation.xlsx")
dir_plots       <- "plots"

if (!dir.exists(dir_plots)) dir.create(dir_plots, recursive = TRUE)

# 2. Load Data & Compute DE Counts ---------------------------------------------

bacteria_results <- read_csv(file_de_results, show_col_types = FALSE)

# Summary counts of DE genes per organism for plot annotation
summary_counts <- bacteria_results %>%
  filter(diffexpressed %in% c("+fiber", "-fiber")) %>%
  group_by(bacteria, diffexpressed) %>%
  summarize(count = n(), .groups = "drop")

# 3. Cayman Substrate Mapping --------------------------------------------------

substrate_mapping <- read_excel(file_cayman, sheet = 1) %>%
  select(Subfamily, Glycan_annotation, `Activity description`) %>%
  mutate(
    substrate = case_when(
      str_detect(Glycan_annotation, "xylan") & str_detect(Glycan_annotation, "mucin_glycan") ~ "both",
      str_detect(Glycan_annotation, "mucin_glycan") ~ "mucin",
      str_detect(Glycan_annotation, "xylan") ~ "xylan",
      TRUE ~ "other"
    )
  ) %>%
  distinct(Subfamily, .keep_all = TRUE)

# 4. Integrate Annotations & Assign Categories ---------------------------------

bacteria_cz <- bacteria_results %>%
  dplyr::rename(kegg_id = kegg_consensus) %>%
  left_join(substrate_mapping, by = c("cayman_family" = "Subfamily")) %>%
  mutate(
    category = case_when(
      exists("galactose_metabolism_kos") & kegg_id %in% galactose_metabolism_kos ~ "galactose_metabolism",
      substrate == "xylan" ~ "xylan_degradation",
      substrate == "mucin" ~ "mucin_degradation",
      substrate == "both" ~ "xylan_and_mucin_degradation",
      exists("xylose_transport_kos") & kegg_id %in% xylose_transport_kos ~ "xylose_transport",
      exists("xylose_metabolism_kos") & kegg_id %in% xylose_metabolism_kos ~ "xylose_metabolism",
      exists("arabinose_transport_kos") & kegg_id %in% arabinose_transport_kos ~ "arabinose_transport",
      exists("arabinose_metabolism_kos") & kegg_id %in% arabinose_metabolism_kos ~ "arabinose_metabolism",
      exists("multiple_sugar_transport_kos") & kegg_id %in% multiple_sugar_transport_kos ~ "multiple_sugar_transport",
      exists("inositol_pathway_kos") & kegg_id %in% inositol_pathway_kos ~ "inositol_pathway",
      exists("simple_sugar_transport_kos") & kegg_id %in% simple_sugar_transport_kos ~ "simple_sugar_transport",
      exists("pts_system_kos") & kegg_id %in% pts_system_kos ~ "pts_system",
      exists("mannose_pts_kos") & kegg_id %in% mannose_pts_kos ~ "pts_system",
      !is.na(Glycan_annotation) & !str_detect(Glycan_annotation, "mucin|xylan") ~ Glycan_annotation,
      TRUE ~ "other"
    )
  ) %>%
  select(-any_of("delabel"))

# 5. Extract Gene Labels for Specific Categories --------------------------------

target_categories <- c(
  "xylan_degradation", "arabinose_metabolism", "xylose_metabolism",
  "galactose_metabolism", "pts_system", "abc_sugar_transporter", "inositol_pathway"
)

gene_labels <- bacteria_cz %>%
  filter(diffexpressed != "NO") %>%
  filter(!is.na(gene) | !is.na(cayman_family) | !is.na(dbcan_family)) %>%
  filter(category %in% target_categories) %>%
  mutate(
    delabel = case_when(
      is.na(cayman_family) & is.na(gene) & exists("master_kegg_dict") & kegg_id %in% names(master_kegg_dict) ~ unname(master_kegg_dict[kegg_id]),
      is.na(cayman_family) ~ gene,
      TRUE ~ cayman_family
    )
  ) %>%
  select(locus_tag, bacteria, delabel)

# Merge back labels & remove duplicates per organism
bacteria_cz <- bacteria_cz %>%
  left_join(gene_labels, by = c("locus_tag", "bacteria")) %>%
  mutate(delabel = replace_na(delabel, "")) %>%
  group_by(bacteria, delabel) %>%
  mutate(delabel = if_else(duplicated(delabel) & delabel != "", "", delabel)) %>%
  ungroup()

# 6. Volcano Plot Generation ---------------------------------------------------

color_palette <- c("NO" = "grey70", "+fiber" = "steelblue3", "-fiber" = "slateblue4")

ggplot(bacteria_cz, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = diffexpressed, fill = diffexpressed), alpha = 0.5, size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # DE Gene Count Annotations per Panel
  geom_text(
    data = summary_counts,
    aes(
      x = if_else(diffexpressed == "+fiber", Inf, -Inf),
      y = Inf,
      label = count,
      color = diffexpressed,
      hjust = if_else(diffexpressed == "+fiber", 1.2, -0.2)
    ),
    vjust = 1.3,
    size = 3.8,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  # Gene Labels
  geom_text_repel(
    aes(label = delabel, color = diffexpressed),
    size = 2.8,
    max.overlaps = 10000,
    box.padding = 0.3,
    point.padding = 0.2,
    show.legend = FALSE
    
  ) +
  
  # Styling & Scales
  facet_wrap(~bacteria, ncol = 2) +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  labs(
    x = expression(Log[2]~"Fold Change"),
    y = expression(-Log[10]~"Adjusted "*italic(P)*"-value")
  ) +
  theme_pubclean() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic", size = 10),
    axis.title = element_text(size = 11, face = "bold")
  )

# Save Plot
ggsave(
  filename = file.path(dir_plots, "volcanoplot_poster_genelabels_xylan_all_labels.svg"),
  #plot = p_volcano,
  width = 8,
  height = 11
)
