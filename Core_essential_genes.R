# Description: Analyze 2017_SBC5 screen to verify if selected genes from MuscKA are essential.
# Goal is to profile the gene hits of muscka and DEG and match their essentiality to TKO database.
# Date created: 2025.09.12
# Author: Bell Wu

# 1.0 setwd and load dfs ----
setwd("~/R_programming/R_MuSCK_Library/CRISPR-Screen BL1 data/Essential genes/")
ess_genes_T1 = read.table("T0_v_T25_Essential_Genes.gene_summary.txt", header = TRUE)
ess_genes_T2 = read.table("T0_v_T35_Essential_Genes.gene_summary.txt", header = TRUE)
core_genes = read.table("core-essential-genes-sym_HGNCID.txt", header = TRUE)

# 2.0 identify essential genes: ------
gene_list = c("TAF7", "DDX23", "HNRNPK", "TICRR", "IARS", "AURKB", "PE01", 
              "SETD1B", "CCT3", "CDC25B", "BIRC5", "NRK", "MRG15", "POP1",
              "MED23")
gene_list[!gene_list %in% core_genes$AAMP]
