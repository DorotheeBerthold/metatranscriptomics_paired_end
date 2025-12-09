library(ggpubr)
library(genbankr)
library(dplyr)

# Import shrunken results from DESEq2 pipeline
results <- read.csv("results/DB092_deseq2_results_shrunken.csv", row.names = 1)
results$genes <- rownames(results)

# separate genes column into two based on .
bacteria_results <- results %>% 
  separate(genes, into = c("gene", "bacteria"), sep = "\\.") %>% 
  filter(!is.na(padj)) %>%
  mutate(diffexpressed = recode(diffexpressed, 
                                "NO" = "", 
                                "UP" = "+fiber", 
                                "DOWN" = "-fiber")) 

bacteria_results <- bacteria_results %>%
  mutate(bacteria = recode(bacteria,
                            
                           "I48" = "Bacteroides caecimuris I48",
                           "YL58" = "Blautia coccoides YL58",
                           "YL32" = "Enterocloster clostridioformis YL32",
                           "YL45" = "Turicimonas muris YL45",
                           "I46" = "Clostridium innocuum I46",
                           "YL44" = "Akkermansia muciniphila YL44",
                           "YL31" = "Flavonifractor plautii YL31",
                           "YL2" = "Bifidobacterium animalis YL2"
  ))


summary_counts <- bacteria_results %>%
  group_by(bacteria, diffexpressed) %>%
  summarise(count = n()) %>%
  ungroup() %>% 
  filter(diffexpressed != "")


# make volcano plot for bacteria
# add line ad 0

ggplot(bacteria_results, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
  geom_point(alpha = 0.5) +
  theme_pubclean() +
  facet_wrap(~bacteria) +
  geom_text(
    data = summary_counts,
    aes(
      x = case_when(
        diffexpressed == "+fiber" ~ Inf,
        diffexpressed == "-fiber" ~ -Inf
      ),
      y = Inf,
      label = count,
      color = diffexpressed
    ),
    hjust = case_when(
      summary_counts$diffexpressed == "+fiber" ~ 1.2,
      summary_counts$diffexpressed == "-fiber" ~ -0.2
    ),
    vjust = 1.2,
    size = 4,
    fontface = "bold"  
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("no" = "grey", "+fiber" = "steelblue3", "-fiber" = "slateblue4")) +
  labs(title = "Differential Expression Volcano Plot per Bacteria",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.title = element_blank(), legend.position = "right")
ggsave("plots/volcano_plots_per_bacteria_pseudocount.png", width = 10, height = 8)
