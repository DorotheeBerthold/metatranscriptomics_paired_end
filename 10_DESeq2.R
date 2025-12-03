{
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(apeglm)
}

# Read in raw data
# Genes should be rownames, samples colnames
genes_wide <- read.csv("results/gene_counts_wide_metadata.csv", row.names = 1)
genes_wide <- genes_wide %>%
  select(-c(media, replicate))

genes_t <- as.data.frame(t(genes_wide))


# Create treatment vector
treatment <- factor(c(rep("no_fiber",4), rep("fiber",4)), levels=c("no_fiber", "fiber"))

# Make "annotation" frame, where rownames are sample names & treatment is the different samples
colData <- data.frame(treatment, row.names = colnames(genes_t))

# Create dds project
dds <- DESeqDataSetFromMatrix(
  countData = genes_t, colData = colData, 
  design = ~ treatment)
dim(dds)

# Filter low count genes
idx <- rowSums(counts(dds, normalized=FALSE) >= 1) >= 3
dds.f <- dds[idx, ]
dim(dds.f)

# Estimation of dispersion
dds.f <- DESeq(dds.f)

# Perform PCA to see if there are any outliers
vsd <- varianceStabilizingTransformation(dds.f, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("treatment"))
pcaData + geom_label(aes(x=PC1,y=PC2,label=treatment))

# Downstream results
res <- results(dds.f)
summary( res )

# Plot estimate of dispersions
plotDispEsts(dds.f)

# Look at coefficients of model
head(coef(dds.f))

# MA plot
# Blue dots are significant logFC
res.lfc <- lfcShrink(dds.f, coef=2, res=res)
DESeq2::plotMA(res.lfc)

# Set tresholds
FDRthreshold = 0.01
logFCthreshold = 1.0
# add a column of NAs
res.lfc$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.01, set as "UP" 
res.lfc$diffexpressed[res.lfc$log2FoldChange > logFCthreshold & res.lfc$padj < FDRthreshold] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.01, set as "DOWN"
res.lfc$diffexpressed[res.lfc$log2FoldChange < -logFCthreshold & res.lfc$padj < FDRthreshold] <- "DOWN"


# Make volcanoplot
ggplot( data = data.frame( res.lfc ) , aes( x=log2FoldChange , y = -log10(padj) , col =diffexpressed ) ) + 
  geom_point() + 
  geom_vline(xintercept=c(-logFCthreshold, logFCthreshold), col="red") +
  geom_hline(yintercept=-log10(FDRthreshold), col="red") +
  scale_color_manual(values=c("blue", "grey", "red"))

# Make results table
table(res.lfc$diffexpressed)

# Generate heatmap
vsd.counts <- assay(vsd)

topVarGenes <- head(order(rowVars(vsd.counts), decreasing = TRUE), 20)
mat  <- vsd.counts[ topVarGenes, ] #scaled counts of the top genes
mat  <- mat - rowMeans(mat)  # centering
pheatmap(mat, show_rownames = T, show_colnames = T)

# write csv
write.csv( res ,'results/DB092_deseq2_results.csv' )
