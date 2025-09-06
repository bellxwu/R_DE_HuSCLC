## Description: R script for cleaning CDX file taken from h4h. Filtering protein coding, adding gene IDs and removing duplicates.
## Author: Bell Wu
## Date created: 2025.08.26

library(tidyverse)
library(dplyr)
library(biomaRt)

setwd(dir = "~/R_programming/R_DE_HuSCLC/")

# 1.0 load df and clean data-----
CDX_df <- read.table("CDX_counts.matrix.txt", header = TRUE)
# clean unneccessary columns
View(CDX_df)
str(CDX_df)
CDX_df <- CDX_df |> 
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length"))
head(CDX_df)
# Clean column names
CDX_cols <- sub("(CDX[0-9]+).*", "\\1", colnames(CDX_df))
CDX_cols <- ave(seq_along(CDX_cols), CDX_cols, FUN = function(i) paste0(CDX_cols[i], ".", seq_along(i)))
colnames(CDX_df) = CDX_cols
head(CDX_df)
colnames(CDX_df)[1] <- "ensembl_gene_id"

# dplyr alternative is to convert array to dataframe and the

# 2.0 Identify and filter protein coding genes -----
# retrieve protein coding genes from ensemb
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
                              mart = ensembl)
# next want to clean up this data. Removes ensembl IDs with no HGNC symbols.
dim(protein_coding_genes)
protein_coding_genes <- protein_coding_genes %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "")
dim(protein_coding_genes)

# remove unneeded columns
protein_coding_genes <- protein_coding_genes %>%
  dplyr::select(-external_gene_name, -gene_biotype, -description)
head(protein_coding_genes) # result is df with list of gene_id and hgnc symbol

# check initial dimensions
dim(protein_coding_genes) # check initial dimensions 
head(protein_coding_genes)
dim(CDX_df)
head(CDX_df)

# select only protein coding genes from rows in CDX_df
CDX_protein_coding <- CDX_df[CDX_df$ensembl_gene_id %in% protein_coding_genes$ensembl_gene_id, ] # select out the rows that match the protein coding
dim(CDX_protein_coding) # check dimensions to make sure genes present are at reasonable amt
head(CDX_protein_coding)

# add the gene ID to the database
CDX_protein_coding <- inner_join(CDX_protein_coding, protein_coding_genes, by = "ensembl_gene_id")

# determine if there are any duplicates within HGNC column
dup_logical <- duplicated(CDX_protein_coding$hgnc_symbol) # create a logical vector storing all the duplicated values
duplicated_genes <- CDX_protein_coding$hgnc_symbol[dup_logical] # to find the gene that is duplicated
dup_genes <- CDX_protein_coding |> # subset rows with duplicated gene names
  filter(hgnc_symbol == duplicated_genes)

# remove the duplicated row with no gene counts
CDX_protein_coding <- CDX_protein_coding |> # subset rows with duplicated gene names
  filter(!ensembl_gene_id == "ENSG00000258724")

# write .csv for all cdx read counts
write.csv(CDX_protein_coding, file = "CDX_raw_read_counts.csv")
