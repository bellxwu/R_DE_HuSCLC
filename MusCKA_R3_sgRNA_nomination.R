# Description: Identify most efficient sgRNAs from muscka screen. Goal is to nominate 
# a few sgRNA sequences to use for future KO.
# Date created: 2025.09.09
# Author: Bell Wu
# Changelog:
# 2025.09.09 - nominate sgRNAs for top 5 genes. 
# 2025.09.28- included sgRNAs for new genes based on prevalence

library(tidyverse)

# 1.0 prepare dfs and workspace ------------------------------------------------
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/QC_initialCounts/")

# load data files:
counts_summary <- read.table("MusCK3_r3_ver2.countsummary.txt", header = TRUE)
counts <- read.table("MusCK3_r3_ver2.count.txt", header = TRUE)
counts_norm <- read.table("MusCK3_r3_ver2.count_normalized.txt", header = TRUE)

head(counts_summary)
head(counts)
head(counts_norm)

# 2.0 select out candidate sgRNAs ----------------------------------------------
# raw counts
raw_sgrna = counts |> 
  filter(grepl("GENES", sgRNA)) |> 
  filter(!grepl("GENES", sgRNA)) |> 
  dplyr::arrange(sgRNA)
# normalized
cand_sgrnas = counts_norm |> 
  filter(grepl("GENES", sgRNA)) |> 
  filter(!grepl("GENES", sgRNA)) |> 
  dplyr::arrange(sgRNA)
# select fold change
FC = cand_sgrnas |> 
  mutate(B6_NRG_IgG = log2((B6_IgG+1) / (NRG_IgG+1)),
         B6_NRG_aCD8 = log2((B6_aCD8+1) / (NRG_IgG+1)))
# filter for FC < 1
smal_FC = FC |> 
  filter(B6_NRG_IgG < 1)
top_grna = smal_FC |> 
  group_by(Gene) |> 
  dplyr::arrange(B6_NRG_IgG, .by_group = TRUE)
# take the first two
keep = ((seq_along(top_grna$sgRNA) - 1) %% 5) < 2
select_grnas = top_grna$sgRNA[keep]

# 3.0 load and get sequences ---------------------------------------------------
setwd("~/R_programming/R_MuSCK_Library/")
grna_seq = read.table("liu_lab_musck_library_a_grna_sequences.txt", header = TRUE)
head(grna_seq)
# select out sequences
select_seq = grna_seq |> 
  filter(sgID %in% select_grnas)
# write.csv file
setwd("~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/S4_2025 MusCKA_R3 analysis/")
write.csv(select_seq, file = "Candidate_sgRNAs.csv")

# 4.0 sgRNAs from new candidate genes (based on prevelance) --------------------
new_genes = c("GENES") # OMG!
# filter out genes
cand_sgrnas = counts_norm |> 
  filter(grepl("GENES", sgRNA)) |> 
  dplyr::arrange(sgRNA)
# select fold change
FC = cand_sgrnas |> 
  mutate(B6_NRG_IgG = log2((B6_IgG+1) / (NRG_IgG+1)),
         B6_NRG_aCD8 = log2((B6_aCD8+1) / (NRG_IgG+1)))
# filter for FC < 1, for negatively selected
smal_FC = FC |> 
  filter(B6_NRG_IgG < 1)
# select out top sgRNAs
top_n_grna = smal_FC |> 
  group_by(Gene) |> 
  dplyr::arrange(B6_NRG_IgG, .by_group = TRUE) |> 
  filter(grepl("GENES", sgRNA))
# filter for FC > 1, for positively selected
large_FC = FC |> 
  filter(B6_NRG_IgG > 1)
# select out top sgRNAs
top_p_grna = large_FC |> 
  group_by(Gene) |> 
  dplyr::arrange(B6_NRG_IgG, .by_group = TRUE)
# take the first two (from neg)
keep = ((seq_along(top_n_grna$sgRNA) - 1) %% 5) < 2
select_n_grnas = top_n_grna$sgRNA[keep]
# take the first two (from pos)
keep = ((seq_along(top_p_grna$sgRNA) - 1) %% 5) < 2
select_p_grnas = top_p_grna$sgRNA[keep]
# combine arrays:
select_grnas = c(select_n_grnas, select_p_grnas)

# 4.1  load and get sequences ---------------------------------------------------
setwd("~/R_programming/R_MuSCK_Library/")
grna_seq = read.table("liu_lab_musck_library_a_grna_sequences.txt", header = TRUE)
head(grna_seq)
# select out sequences
select_seq = grna_seq |> 
  filter(sgID %in% select_grnas)
# write.csv file
setwd("~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/S4_2025 MusCKA_R3 analysis/")
write.csv(select_seq, file = "New_Candidate_sgRNAs.csv")


