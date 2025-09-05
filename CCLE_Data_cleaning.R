## Description: Compiling the datasets from CCLE and editing to include the cell-lines of interest
## Author: Bell Wu
## Last updated: 2025.02.06

library(tidyverse)
library(dplyr)
library(biomaRt)

# 1.1: Used for direct data (.gct file) from CCLE. ----
# Set working directory:
getwd()
setwd("~/R_programming/R_DE_HuSCLC/")

# read .gz file taken directly from CCLE
CCLE_counts <- read.delim(gzfile("CCLE_RNAseq_genes_counts_20180929.gct.gz"),
                          skip = 2, header = TRUE, sep = "\t")
CCLE_info <- readLines(gzfile("CCLE_RNAseq_genes_counts_20180929.gct.gz"), n = 2) # to keep the header info of CCLE 
subset_CCLE <- CCLE_counts[1:5, 1:5] # take a small sample of the data

# dataset has duplicate names. Different version of gene is found. 
# Worthwhile to first select the cell-lines of interest first before anything else.

selected_line <- CCLE_counts |> 
  dplyr::select(matches("Name|Description|H841|SHP77|H82|H446|SBC5|H1694|H1048|H526"))
head(selected_line)

# #install biomaRt
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")

# use biomaRt to retrieve protein coding genes
ensembl <- useEnsembl(biomart = "genes")
searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# filtering out only ensembl IDs to be protein coding genes
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
searchAttributes(mart = ensembl, pattern = "biotype")
# use the getBM() function to query for only protein coding genes
protein_coding_genes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name", "description", "gene_biotype"),
      filters = "biotype", # select the biotype 
      values = "protein_coding", # select only genes with a protein_coding biotype
      mart = ensembl
      ) # identify protein coding genes

# next want to clean up this data. Removes ensembl IDs with no HGNC symbols.
dim(protein_coding_genes)
protein_coding_genes <- protein_coding_genes %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "")
dim(protein_coding_genes)

# remove the unnecessary selected columns
protein_coding_genes <- protein_coding_genes %>%
  dplyr::select(-external_gene_name, -gene_biotype, -description)
head(protein_coding_genes)

# now that I have an ensembl gene list with only protein coding genes, want to then match with the CCLE database to filter out
dim(protein_coding_genes) # check initial dimensions 
head(protein_coding_genes)
dim(selected_line)
head(selected_line)

# need to separate out the versions as CCLE ensembl ID likely from an outdated version. Remove the version
selected_line_sep <- separate(selected_line, Name, into = c("ensembl", "version"),
                                            sep = "\\.")
selected_line_sep <- selected_line_sep %>%
  dplyr::select(-version) # remove version column

# alternatively can just remove the version using dplyr
selected_line_no_ver <- selected_line |> 
  mutate(Name = sub("\\..*", "", Name))
head(selected_line_no_ver)

# select for the protein coding genes only 
head(protein_coding_genes)
SCLC_protein_coding <- selected_line_no_ver[selected_line_no_ver$Name %in% protein_coding_genes$ensembl_gene_id, ]
dim(SCLC_protein_coding) # check dimensions to make sure genes present are at reasonable amt
head(SCLC_protein_coding)

# determine if there are any duplicates within HGNC column
dup_logical <- duplicated(SCLC_protein_coding$Description) # create a logical vector storing all the duplicated values
duplicated_genes <- SCLC_protein_coding$Description[dup_logical] # to find the gene that is duplicated
duplicated_genes <- SCLC_protein_coding[dup_logical, ] # to subset the entire row that is duplicated

# subset out the rows that are duplicated by gene name
duplicates <- SCLC_protein_coding %>%  
  dplyr::filter(Description == "DEFB130") # selects out all rows with duplicate gene name "DEFB130"

# since values are 0, can just remove this row
dim(SCLC_protein_coding)
SCLC_protein_coding <- SCLC_protein_coding %>%
  dplyr::filter(!Description == "DEFB130") # use the filter function with ! operator to remove all rows that matches this name
dim(SCLC_protein_coding)

# save as .csv
head(SCLC_protein_coding)
write.csv(SCLC_protein_coding, file = "Protein-coding_SCLC_cell_lines.csv", row.names = FALSE)


# 1.2: Used for CCLE file (.csv) where gene names are separated with .. -----
# read csv files
SCLC_CCLE <- read.csv("SCLC_CCLE.csv")
subset_SCLC <- head(SCLC_CCLE) # observe layout

library(tidyr) 
#code with smaller subset
subset_SCLC <- separate(subset_SCLC, X, into = c("Gene", "GeneID"), 
                        sep = "\\.\\.")
row.names(subset_SCLC) = subset_SCLC$Gene
subset_SCLC <- subset_SCLC[ , c(-1,-2)]

# apply code to larger subset
SCLC_CCLE <- separate(SCLC_CCLE, X, into = c("Gene", "GeneID"), 
                        sep = "\\.\\.")
# NAs are detected, will need to identify them
incomplete_rows <- SCLC_CCLE[!complete.cases(SCLC_CCLE), ] # identify which row has the NA
dim(SCLC_CCLE) # check dimensions
SCLC_CCLE <- SCLC_CCLE[complete.cases(SCLC_CCLE), ]
dim(SCLC_CCLE) # check dimensions again to see if it worked
# Great success! - Borat

row.names(SCLC_CCLE) = SCLC_CCLE$Gene
SCLC_CCLE <- SCLC_CCLE[ , c(-1,-2)]

# select cell-lines of interest
selected_lines <- SCLC_CCLE %>%
  select(contains("H841"), contains("SHP77"), contains("H82"), contains("H446"), 
         contains("SBC5"), contains("H1694")
         )

#export as .csv
write.csv(selected_lines, file = "SCLC_cell_lines.csv", row.names = TRUE)

# after running this script, will have a list of cell-lines of interest




