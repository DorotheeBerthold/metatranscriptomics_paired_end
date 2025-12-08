# DESEq2 normalisation ===================================================================
# Code according to Klingenberg and Meinicke 2017:
# "How to normalize metatranscriptomic count data for differential expression analysis"
# ========================================================================================



# Libraries ------------------------------------------------------------------------------

library(rgdal)      
library(ggplot2)    
library(ggrepel)    
library(data.table)
library(tidyverse)
library(DESeq2)
library(abind)
library(pheatmap)

# Save pheatmap function -----------------------------------------------------------------

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Data import ----------------------------------------------------------------------------

abundances_TX <- read.csv("results/abundances_TX_sorted.csv", row.names = 1)
# Remove Muribaculum & Lactobacillus columns
abundances_TX <- abundances_TX %>%
  select(-c(Muribaculum.intestinale.YL27, Lactobacillus.reuteri.I49))
genes_wide <- read.csv("results/gene_counts_wide_metadata.csv", row.names = 1)
genes_wide <- genes_wide %>%
  select(-c(media, replicate))
genes_t <- as.data.frame(t(genes_wide))

# remove genes from YL27 & I49
genes_t <- genes_t %>% 
  filter(!grepl("YL27", rownames(genes_t)),
         !grepl("I49", rownames(genes_t)))

# Split big gene matrix by organism -------------------------------------------------------

# Add 1 to every count to avoid zeros
genes_t <- genes_t + 1

# extract organism ID from row names
genes_t$org <- sub(".*\\.", "", rownames(genes_t)) 

# Filter out rows where rowSums are <300
#genes_t <- genes_t[rowSums(genes_t[, 1:ncol(genes_t)-1]) >= 300, ]



# list of organism IDs present in count table
org_list <- unique(genes_t$org)

# split into count matrices per organism
all_org_counts <- lapply(org_list, function(org){
  mat <- genes_t[genes_t$org == org, 1:ncol(genes_t)-1]  # remove org col
  rownames(mat) <- rownames(genes_t)[genes_t$org == org]
  return(as.matrix(mat))
})

names(all_org_counts) <- org_list

# Extract dimnames from first dimension from all_org_padded as list
Xlist <- lapply(all_org_counts, function(mat) mat[,1 ])

# pad to make them all same length
pad_to <- max(sapply(all_org_counts, nrow))

# Create vector to store number of padded rows per organism
padded_rows <- integer(length(all_org_counts)) 

all_org_counts_padded <- lapply(seq_along(all_org_counts), function(i){
  mat <- all_org_counts[[i]]
  n_pad <- pad_to - nrow(mat)
  padded_rows[i] <<- n_pad
  if (n_pad > 0) {
    pad <- matrix(0L, nrow = n_pad, ncol = ncol(mat))
    mat <- rbind(mat, pad)
  }
  mat
})

names(all_org_counts_padded) <- org_list

# convert abundance table to matrix
Sampling.Rate.Mat <- as.matrix(abundances_TX)

# extract organism IDs from abundance table
colnames(Sampling.Rate.Mat) <- sub(".*\\.", "", colnames(Sampling.Rate.Mat))

# ensure the order of the organims match
org_order <- names(all_org_counts_padded)
Sampling.Rate.Mat <- Sampling.Rate.Mat[, org_order]


# Function to build 3D array using abind -----------------------------------------------------

abind.matrix <- function(Xarray=NULL,XMat)
{
  require(abind)
  #Xmat = raw count data for one species
  #Xmat rows = features, coloums = samples
  #Xarray = array containing data for already added species
  Xarray <- abind(Xarray,XMat, along = 3)
  return(Xarray)
}

# Build 3D array of raw counts ----------------------------------------------------------

Xarray <- NULL
for (org in org_order) {
  XMat <- all_org_counts_padded[[org]]
  Xarray <- abind.matrix(Xarray, XMat)
}
dim(Xarray)  # check dimensions (genes x samples x organisms)


# Create condition vector ----------------------------------------------------------------

cond.vec <- factor(c(rep("no_fiber",4), rep("fiber",4)))

# Factor so that "UP" is higher in fiber condition
cond.vec <- factor(cond.vec, levels = c("no_fiber", "fiber"))

# DESeq2 normalisation function ----------------------------------------------------------

DESeq2.norm.mat <- function(Xmat, cond.vec) {
  # Round counts and ensure integer
  Xmat <- round(Xmat)
  storage.mode(Xmat) <- 'integer'
  
  # Build DESeq2 dataset
  colData <- data.frame(condition = factor(cond.vec, levels = unique(cond.vec)))
  dds <- DESeqDataSetFromMatrix(countData = Xmat, colData = colData, design = ~condition)

  # Normalize counts by size factors
  # Zero-safe normalization
  dds <- estimateSizeFactors(dds, type = "poscounts")
  
  # Get normalized counts
  YMat <- counts(dds, normalized = TRUE)
  
  return(YMat)
}

# DESeq2 taxon-specific normalisation function ----------------------------------------------------

DESeq2.tax.specific <- function(Xarray, Xlist, cond.vec) {
  require(DESeq2)
  
  nOrg <- dim(Xarray)[3]
  nSamples <- dim(Xarray)[2]
  
  # Check input
  if (length(cond.vec) != nSamples) {
    stop("Length of cond.vec does not match number of samples")
  }
  
  # Build a list of normalized matrices for each organism
  norm.list <- vector("list", nOrg)
  for (i in seq_len(nOrg)) {
    mat <- Xarray[,,i]
    
    n_pad <- padded_rows[i]
    if (n_pad > 0) {
      mat <- mat[1:(nrow(mat) - n_pad), , drop = FALSE]
    } else {
      # No padding, keep matrix as-is
      mat <- mat
    }
  
  cat("\nOrganism:", i,
      "\nRows after filtering:", nrow(mat),
      "\nGene names available:", length(Xlist[[i]]),
      "\nGene name labels:", length(names(Xlist[[i]])),
      "\n")
  
  # ---- DEBUG CHECKS ----
  if (nrow(mat) == 0) stop(paste("Organism", i, "has zero non-zero rows"))
  if (length(names(Xlist[[i]])) == 0) stop(paste("Xlist[[", i, "]] has NO names"))
  if (nrow(mat) > length(names(Xlist[[i]]))) {
    stop(paste("Organism", i, "has MORE non-zero rows (", nrow(mat),
               ") than gene names (", length(names(Xlist[[i]])), ")"))
  }
    
    # Assign gene names safely
    rownames(mat) <- names(Xlist[[i]])[seq_len(nrow(mat))]
    
    mat <- mat[rowSums(mat) > 0, ]  # further filter zero-sum rows
    
    # Normalize counts
    mat <- DESeq2.norm.mat(mat, cond.vec = cond.vec)
    
    norm.list[[i]] <- mat
    
    # save data frame
    cat("Normalized matrix dimensions for organism", i, ":", dim(mat), "\n")
    
}
  # ---- REPAD ALL NORMALIZED MATRICES BEFORE SUMMING ----
  
  all_genes <- unique(unlist(lapply(norm.list, rownames)))
  
  norm.list.full <- lapply(norm.list, function(mat) {
    missing <- setdiff(all_genes, rownames(mat))
    if (length(missing) > 0) {
      pad <- matrix(0, nrow = length(missing), ncol = ncol(mat))
      rownames(pad) <- missing
      mat <- rbind(mat, pad)
    }
    mat[all_genes, , drop = FALSE]
  })
 
  # Combine all normalized matrices by summing
  Scaled.Mat <- Reduce("+", norm.list.full)
  
  # Ensure valid and unique rownames
  if(any(duplicated(rownames(Scaled.Mat)))) {
    rownames(Scaled.Mat) <- make.unique(rownames(Scaled.Mat))
  }
  
  # Ensure column names are valid
  if(is.null(colnames(Scaled.Mat))) {
    colnames(Scaled.Mat) <- paste0("sample_", seq_len(nSamples))
  }
  
  # Build colData
  colData <- data.frame(condition = factor(cond.vec, levels = unique(cond.vec)))
  rownames(colData) <- colnames(Scaled.Mat)
  
  # Round counts and convert to integer for DESeq2
  Scaled.Mat <- round(Scaled.Mat)
  storage.mode(Scaled.Mat) <- "integer"
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = Scaled.Mat,
    colData = colData,
    design = ~condition
  )
  
  # Zero-safe size factor estimation
  dds <- estimateSizeFactors(dds, type = "poscounts")
  
  return(list(norm.list = norm.list, dds = dds, Scaled.Mat = Scaled.Mat))
}


# DESeq2 taxon-specific normalisation ----------------------------------------------------

dds <- DESeq2.tax.specific(
  Xarray     = Xarray,
  cond.vec   = cond.vec,
  Xlist = Xlist

)

# extract norm list from dds & dds
norm.list <- dds$norm.list
names(norm.list) <- org_list
#saveRDS(norm.list, file = "results/DESeq2_taxon_specific_normalized_counts.rds")

# extract normalised matrix
scaled_mat <- dds$Scaled.Mat
#write.csv(scaled_mat, file = "results/DESeq2_taxon_specific_scaled_matrix.csv")

# extract DESEq2 object
dds <- dds$dds

dim(dds)

# Filter low count genes
idx <- rowSums(counts(dds, normalized=FALSE) >= 2) >= 3
dds.f <- dds[idx, ]
dim(dds.f)

# Estimation of dispersion
dds.f <- DESeq(dds.f)

# Perform PCA to see if there are any outliers
vsd <- varianceStabilizingTransformation(dds.f, blind=TRUE )
pcaData <- plotPCA(vsd, intgroup=c("condition"))
pcaData + geom_label(aes(x=PC1,y=PC2,label=condition))

# Downstream results
res <- results(dds.f)
summary( res )

# Plot estimate of dispersions
plotDispEsts(dds.f)
ggsave("plots/dispersion_dds.png")

# Look at coefficients of model
head(coef(dds.f))

# MA plot
# Blue dots are significant logFC
res.lfc <- lfcShrink(dds.f, coef=2, res=res)
DESeq2::plotMA(res.lfc)

res_log2FC <- res.lfc$log2FoldChange


# Set tresholds
FDRthreshold = 0.01
logFCthreshold = 1.0
# add a column of NAs
res.lfc$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.01, set as "UP" 
res.lfc$diffexpressed[res.lfc$log2FoldChange > logFCthreshold & res.lfc$padj < FDRthreshold] <- "UP"
# if log2Foldchange < 1 and pvalue < 0.01, set as "DOWN"
res.lfc$diffexpressed[res.lfc$log2FoldChange < -logFCthreshold & res.lfc$padj < FDRthreshold] <- "DOWN"

# Filter NA padj for plotting
res.plot <- res.lfc[!is.na(res.lfc$padj), ]
# Make volcanoplot
ggplot( data = data.frame( res.plot ) , aes( x=log2FoldChange , y = -log10(padj) , col =diffexpressed ) ) + 
  geom_point() + 
  geom_vline(xintercept=c(-logFCthreshold, logFCthreshold), col="red") +
  geom_hline(yintercept=-log10(FDRthreshold), col="red") +
  scale_color_manual(values=c("blue", "grey", "red"))

# Make results table
table(res.lfc$diffexpressed)

# Generate heatmap
vsd.counts <- assay(vsd)

# Most variable genes heatmap

topVarGenes <- head(order(rowVars(vsd.counts), decreasing = TRUE), 20)
mat  <- vsd.counts[ topVarGenes, ] #scaled counts of the top genes
mat  <- mat - rowMeans(mat)  # centering
mat <- pheatmap(mat, show_rownames = T, show_colnames = T)

save_pheatmap_pdf(mat,
  filename = "plots/top_variable_genes_heatmap.pdf",
  width = 7,
  height = 5
)

# most DE genes heatmap
topDE <- head(order(res$padj), 20)
mat <- vsd.counts[topDE, ]
mat <- pheatmap(mat, show_rownames = T, show_colnames = T)

save_pheatmap_pdf(mat,
                  filename = "plots/top_differential_genes_heatmap.pdf",
                  width = 7,
                  height = 5
)

# write csv
write.csv( res ,'results/DB092_deseq2_results.csv' )
write.csv(as.data.frame(res.lfc), "results/DB092_deseq2_results_shrunken.csv")

