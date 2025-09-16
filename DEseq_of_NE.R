# Description: DE of subtypes (NE-low and NE-high). Includes CDX in this analysis.
# Author: Bell Wu
# Date created: 2025.08.30

library(DESeq2)
library(tidyverse)

setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")

# 1.0 Load files and prep dfs -------
total_counts <- read.csv("CDX_CCLE_raw_counts.csv")
head(total_counts)

# subset out only NE-low, NE-high
NE_low <- total_counts |> 
  dplyr::select(matches("CDX12|H1048|SBC5|H1048|H526"))
head(NE_low)
rownames(NE_low) = total_counts$X

NE_high <- total_counts |> 
  dplyr::select(matches("CDX15|CDX9|H1694|H446|H82|SHP77"))
head(NE_high)
rownames(NE_high) = total_counts$X

# 2.0 DEseq for NE-high:
# 2.1 create metadata for dfs -----
colnames(total_counts)
batches <- data.frame(
  Sample = colnames(total_counts),
  batch  = ifelse(grepl("LUNG$", colnames(total_counts)), "CCLE",
                  ifelse(grepl("^CDX",  colnames(total_counts)), "CDX", NA))
)
NE_high_md <- batches |> 
  filter(batches$Sample %in% colnames(NE_high)) |> 
  mutate(condition = ifelse(
    str_detect(Sample, "CDX9|SHP77|H1694"), "Hot", "Cold"
  ))

# 2.2 Create DDS for DESeq2 -----
high_dds <- DESeqDataSetFromMatrix(countData = NE_high,
                                  colData = NE_high_md,
                                  design = ~ batch + condition)
high_dds <- DESeq(high_dds) # identify DESeq 
high_res <- results(high_dds) # show results
str(high_res)
# identify comparisons
resultsNames(high_dds)
high_resLFC <- lfcShrink(high_dds, coef = "condition_Hot_vs_Cold")
view(high_resLFC)

# probing results
resOrdered <- high_res[order(high_res$pvalue),]
sum(high_res$padj < 0.05, na.rm = TRUE) # how many padj < 0.05
high_resSig <- subset(resOrdered, padj < 0.05)

write.csv(as.data.frame(high_resSig), 
          file="NE_high_HotvsCold.csv")

# 2.3 Visualizing results and QC: -----
plotMA(high_resLFC)
plotMA(high_res)

# 3.0 DEseq for NE-low:
# 3.1 create metadata for dfs -----
NE_low_md <- batches |> 
  filter(batches$Sample %in% colnames(NE_low)) |> 
  mutate(condition = ifelse(
    str_detect(Sample, "CDX12|SBC5|H841"), "Hot", "Cold"
  ))
# 3.2 Create DDS for DESeq2 -----
low_dds <- DESeqDataSetFromMatrix(countData = NE_low,
                                   colData = NE_low_md,
                                   design = ~ batch + condition)
low_dds <- DESeq(low_dds) # identify DESeq 
low_res <- results(low_dds) # show results
str(low_res)
# identify comparisons
resultsNames(low_dds)
low_resLFC <- lfcShrink(low_dds, coef = "condition_Hot_vs_Cold")

# probing results
resOrdered <- low_res[order(low_res$pvalue),]
sum(low_res$padj < 0.05, na.rm = TRUE) # how many padj < 0.05
low_resSig <- subset(resOrdered, padj < 0.05)

write.csv(as.data.frame(low_resSig), 
          file="NE_low_HotvsCold.csv")

# 4.0 Comparisons for same genes ----
NE_low_sig <- read.csv("NE_low_HotvsCold.csv")
NE_high_sig <- read.csv("NE_high_HotvsCold.csv")
head(NE_low_sig)
dim(NE_low_sig)
head(NE_high_sig)
dim(NE_high_sig)
# identify same genes
same_genes <- intersect(NE_low_sig$X, NE_high_sig$X)
length(same_genes)
# subset out
same_csv_NE_low <- NE_low_sig |> 
  filter(NE_low_sig$X %in% same_genes)

same_csv_NE_high <- NE_high_sig |> 
  filter(NE_high_sig$X %in% same_genes)

# write csv
write.csv(same_csv_NE_high, "NE_high_sg.csv")
write.csv(same_csv_NE_low, "NE_low_sg.csv")

# 5.0 probing through same genes ----
view(same_csv_NE_high)
view(NE_high_sig)
view(NE_low_sig)

