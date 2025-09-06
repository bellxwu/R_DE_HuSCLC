## Description: Checking for batch effect in CDX vs CCLE database
## Author: Bell Wu
## Date created: 2025.08.29

library(dplyr)
library(DESeq2)
library(limma)

setwd("~/R_programming/R_DE_HuSCLC/")

# load data files:
CDX <- read.csv("CDX_raw_read_counts.csv")
cell_line <- read.csv("Protein-coding_SCLC_cell_lines.csv")

# filter for only CDX12
head(CDX)
CDX <- CDX |> 
  dplyr::select(c("CDX12.1", "CDX12.2", "hgnc_symbol", "ensembl_gene_id", "CDX15.1", "CDX9.1", "CDX9.2"))
# change column name of ensemble ID
colnames(CDX)[4] <- "Name"
head(cell_line)
# join the dfs
total_counts <- inner_join(CDX, cell_line, by = "Name")
head(total_counts)
# create a .csv
write.csv(total_counts, "CDX_CCLE_raw_counts.csv")

# create batch metadata df
colnames(total_counts)
batches <- data.frame(
  Sample = colnames(total_counts),
  batch  = ifelse(grepl("LUNG$", colnames(total_counts)), "CCLE",
                  ifelse(grepl("^CDX",  colnames(total_counts)), "CDX", NA)),
) |> 
  factor(batch)
# remove all rows that contain NA
batches <- batches |> 
  filter(!is.na(batches$batch))
batches <- batches |> 
  mutate(condition = 
           case_when(str_detect(Sample, "CDX12|H1048|H526|H841|SBC5") ~ "NE_low",
                     str_detect(Sample, "CDX15|CDX9|H1694|H446|H82|SHP77") ~ "NE_high"))

# normalize: ----
# need to remove character values
head(total_counts)
genes <- total_counts$hgnc_symbol
total_counts <- total_counts |> 
  dplyr::select(-c("hgnc_symbol", "Name", "Description"))
# add gene names to row
rownames(total_counts) = genes

# normalize with vst
batch_dds <- DESeqDataSetFromMatrix(countData = total_counts,
                                    colData = batches,
                                    design = ~ batch + condition)
batch_vst <- vst(batch_dds, nsub = 5000)

# plot pca
plotPCA(batch_vst, intgroup = c("batch", "condition"))








