## Description: R script for writing differential analysis of selected SCLC cell-lines
## Author: Bell Wu
## Date Created: 2025.02.06

library(tidyverse)
library(dplyr)
library(biomaRt)
library(DESeq2)

# 1.1: Setting up workspace

# Set working directory:
setwd("~/R_programming/R_DE_HuSCLC/")

# Read .csv of the selected cell-lines from CCLE
selected_lines <- read.csv("Protein-coding_SCLC_cell_lines.csv")
dim(selected_lines)
row.names(selected_lines) = selected_lines$Description
selected_lines <- selected_lines[ , c(-1:-2)]
head(selected_lines)

# need to make a colData containing the metadata and information about the samples
# 1.2 N vs A ----
cell_lines <- selected_lines %>%
  select(contains("H82"), contains("H446"), contains("SHP77"), contains("H1694")) # this will be my input count matrix
head(cell_lines)

sample_names <- colnames(cell_lines)
sample_conditions <- c("Desert", "Desert", "Excluded", "Excluded") # in accordance to the samples themselves
sampleInfo <- data.frame(condition = factor(sample_conditions)) # create a factorized dataframe of the conditions
row.names(sampleInfo) = sample_names
sampleInfo # to check the created data frame

# once have the counts and sample info, can create the dds
dds <- DESeqDataSetFromMatrix(countData = cell_lines,
                              colData = sampleInfo,
                              design = ~ condition)
dds <- DESeq(dds)
DEG_DvE <- results(dds)
head(DEG_DvE)

# can standardize (or scale) the results using LFC
resultsNames(dds)
DEG_DvE_LFCapel <- lfcShrink(dds, coef = "condition_Excluded_vs_Desert", type = "apeglm") 
# shrink using the apeglm method, there are multiple methods specified. Can try them differently to see results. 
# For this purpose will go with the apeglm which is standardly used. 

# ordering the pvalues
DEG_DvE_ordered <- DEG_DvE[order(DEG_DvE$pvalue), ] # w/o LFC shrinkage
DEG_DvE_LFCapel_ordered <- DEG_DvE_LFCapel[order(DEG_DvE_LFCapel$pvalue), ] # w/ LFC shrinkage

# summarizing the list of hits. 
summary(DEG_DvE_ordered)
summary(DEG_DvE_LFCapel_ordered)

sum(DEG_DvE$padj < 0.1, na.rm=TRUE) # identify all genes that have a padj < 0.1


# exporting to a .csv
write.csv(as.data.frame(DEG_DvE_ordered), 
          file = "SCLC_CCLE_Cold-D_vs_Cold-E.csv") # this final .csv contains all genes ordered from most to least significant

# 1.3 I vs P
cell_lines <- selected_lines %>%
  dplyr::select(matches("H1048|H841|SBC5|H526")) # this will be my input count matrix
head(cell_lines)

sample_names <- colnames(cell_lines)
sample_conditions <- c("SCLC-P", "SCLC-P", "Inflamed", "Inflamed") # in accordance to the samples themselves
sampleInfo <- data.frame(condition = factor(sample_conditions)) # create a factorized dataframe of the conditions
row.names(sampleInfo) = sample_names
sampleInfo # to check the created data frame

# once have the counts and sample info, can create the dds
dds <- DESeqDataSetFromMatrix(countData = cell_lines,
                              colData = sampleInfo,
                              design = ~ condition)
sdds <- DESeq(dds)
DEG_DvE <- results(sdds)
head(DEG_DvE)

# can standardize (or scale) the results using LFC
resultsNames(dds)
DEG_DvE_LFCapel <- lfcShrink(dds, coef = "condition_Excluded_vs_Desert", type = "apeglm") 
# shrink using the apeglm method, there are multiple methods specified. Can try them differently to see results. 
# For this purpose will go with the apeglm which is standardly used. 

# ordering the pvalues
DEG_DvE_ordered <- DEG_DvE[order(DEG_DvE$pvalue), ] # w/o LFC shrinkage
DEG_DvE_LFCapel_ordered <- DEG_DvE_LFCapel[order(DEG_DvE_LFCapel$pvalue), ] # w/ LFC shrinkage

# summarizing the list of hits. 
summary(DEG_DvE_ordered)
summary(DEG_DvE_LFCapel_ordered)

sum(DEG_DvE$padj < 0.1, na.rm=TRUE) # identify all genes that have a padj < 0.1


