library(tidyverse)

raw <- readLines("summary_rRNA.txt")

# Extract relevant lines
mRNA      <- raw[grep("mRNA pairs kept", raw)]
rRNA      <- raw[grep("rRNA pairs kept", raw)]
samples <- sub("^Processing\\s+", "", raw[grep("^Processing", raw)])


df <- tibble(
  sample = samples,
  mRNA = as.integer(str_extract(mRNA, "\\d+")),
  rRNA = as.integer(str_extract(rRNA, "\\d+"))
)

# Plot stacked barplot for rRNA removal
df_long <- df %>%
  pivot_longer(cols = c("mRNA", "rRNA"), names_to = "status", values_to = "count")

library(ggplot2)
ggplot(df_long, aes(x = sample, y = count, fill = status)) +
  geom_bar(stat = "identity") +
  labs(title = "rRNA Removal Overview", x = "", y = "Number of Reads", fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/rRNA_removal_overview_update.png", width = 10, height = 6)

# Plot percentage of rRNA vs mRNA
# calculate percentage
df$total <- df$mRNA + df$rRNA

df <- df %>%
  mutate(
    mRNA_perc = round((mRNA / total) * 100, 2),
    rRNA_perc = round((rRNA / total) * 100, 2)
  )

df_perc <- df %>%
  select(sample, mRNA_perc, rRNA_perc) %>%
  pivot_longer(cols = c("mRNA_perc", "rRNA_perc"), names_to = "status", values_to = "percentage")


ggplot(df_perc, aes(x = sample, y = percentage, fill = status)) +
  geom_bar(stat = "identity") +
  labs(title = "rRNA Removal Overview", x = "", y = "Percentage of total reads", fill = "Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("plots/rRNA_removal_percentage_update.png", width = 10, height = 6)
