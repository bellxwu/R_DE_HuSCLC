# Description: Combination analysis from DESeq2 HotvsCold and MusCKA IC vs ID.
# Takes muscka hits with DEGs and finds intersecting genes in DEG list.
# Also filters out muscka hits found in B6-aCD8 vs NRG-IgG.
# Author: Bell Wu
# Date created: 2025.09.03

# Changelog:
# 2025.09.03 - Analysis as described above
# 2025.09.10 - Added t1, t2, t3 analyses. See companion .md "mageck and de data analysis" for rationale and implications of these analyses. 
# 2025.09.12 - Analysis by reframing sign of lfc shown added in 4.0. Instead of "same sign analysis," will now include 
# 2025.09.17 - Created .csv for the gene lists of t1 vs t3; and NE-high + muscka, NE-low + muscka
# opposite sign by matching neg.p of muscka with pos.lfc of DEGs. See companion .md for rationale.

library(tidyverse)

# 1.0 B6 IgG vs NRG IgG, or "t1" ----------
# loss of these genes result in sensitization in IC setting
# load dfs
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")

NE_low_sig <- read.csv("NE_low_HotvsCold.csv")
NE_high_sig <- read.csv("NE_high_HotvsCold.csv")
NE_total <- rbind(NE_high_sig, NE_low_sig)
head(NE_total)
# view(NE_total)

setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
B6_NRG <- read.table("B6_IgG_vs_NSG_IgG_gene_summary.txt", header = TRUE)
head(B6_NRG)
t1_neg_rank100 <- B6_NRG |> 
  filter(neg.rank <= 100)
dim(t1_neg_rank100)

# 1.1 Comparison of DE and muscka
# select genes
head(NE_total)
imm_induce_DE <- NE_total |>
  filter(log2FoldChange > 1)
dim(imm_induce_DE)

head(B6_NRG)
imm_induce_Mus <- B6_NRG |>
  filter(neg.fdr < 0.2)
dim(imm_induce_Mus)

# see if there are any intersecting genes for all three NE low DEs, NE high DEs, and MusCKA
intersect(imm_induce_Mus$id, imm_induce_DE$X)

# see if there are any intersecting genes for either NE_low or NE_high and MusCKA only
common_imm <- intersect(t1_neg_rank100$id, NE_total$X)
DE_genes <- NE_total |> 
  filter(X %in% common_imm)
DE_genes <- DE_genes[order(DE_genes$log2FoldChange, decreasing = TRUE), ]
dim(DE_genes)
mus_DE <- t1_neg_rank100 |> 
  filter(id %in% DE_genes$X)

# intersecting genes in all NE and MusCKA
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
NE_low_sig <- read.csv("NE_high_sg.csv")
NE_high_sig <- read.csv("NE_low_sg.csv")
head(NE_low_sig)
head(NE_high_sig)
intersect(NE_high_sig$X, t1_neg_rank100$id)

# intersecting genes for NE high and MusCKA
length(intersect(NE_high_sig$X, t1_neg_rank100$id))
length(intersect(NE_low_sig$X, t1_neg_rank100$id))

# 2.0 B6-aCD8 vs B6-IgG or "t2": --------
# loss of these genes implies a pro-tumour function of CD8
# alternatively, pos enrichment leads to improved survival in CD8-present environment
setwd("~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/")
aCD8_IgG <- read.table("B6_aCD8_vs_B6_IgG_gene_summary.txt", header = TRUE)
t2_pos_rank100 <- aCD8_IgG |> 
  filter(pos.rank <= 100)
t2_neg_rank100 <- aCD8_IgG |> 
  filter(neg.rank <= 100)

# 2.1 Comparison of muscka t1 vs t2
# see if positive rank for B6-CD8 intersects with neg ranked for B6-IgG
sg <- intersect(t2_pos_rank100$id, t1_neg_rank100$id)
pos_sg <- t2_pos_rank100 |> 
  filter(id %in% sg)
pos_sg <- pos_sg[order(pos_sg$pos.rank) ,]

neg_sg <- t2_neg_rank_100 |> 
  filter(id %in% sg)

# 3.0 B6-aCD8 vs NRG-IgG or "t3" -----
# loss of these genes means cancer becomes sensitized to non-CD8 cells.
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
aCD8_NRG_IgG <- read.table("B6_aCD8_vs_NSG_IgG_gene_summary.txt", header = TRUE)
t3_neg_rank100 <- aCD8_NRG_IgG |> 
  filter(neg.rank <= 100)
t3_neg_0.2 = aCD8_NRG_IgG |> 
  filter(neg.fdr < 0.2)
head(t3_neg_rank100)
t3_pos_rank100 <- aCD8_NRG_IgG |> 
  filter(pos.rank <= 100)

# 3.1 Comparison of muscka t1 vs t3
sg_t3 <- intersect(t3_neg_rank100$id, t1_neg_rank100$id) # these genes are elim regardless of CD8
diff_t1 <- setdiff(t1_neg_rank100$id, t3_neg_rank100$id) # find genes in t1 that are not in t3
length(diff_t1)
# take this from the gene summary table
mus_genes <- t1_neg_rank100 |> 
  filter(id %in% diff_t1)
mus_genes_0.2 <- mus_genes |> 
  filter(neg.fdr <= 0.2)

# 3.2 Compare the filtered genes with the DE list. These should be genes with no NK confounder
sg_mus_DE <- intersect(mus_genes$id, NE_total$X)
length(sg_mus_DE)

mus_genes <- mus_genes |> 
  filter(mus_genes$id %in% sg_mus_DE)
mus_genes <- mus_genes[order(mus_genes$id), ]
DE_genes <- NE_total |> 
  filter(X %in% sg_mus_DE)
DE_genes <- DE_genes[order(DE_genes$X), ]
DE_genes <- DE_genes[order(DE_genes$log2FoldChange, decreasing = TRUE), ]

# 4.0 Finding intersecting genes ------
# take the difference between t1 and t3. 
diff_t1 <- setdiff(t1_neg_rank100$id, t3_neg_rank100$id)
mus_genes <- t1_neg_rank100 |> 
  filter(id %in% diff_t1)
other_mus = t1_neg_rank100 |> 
  filter(!id %in% diff_t1)
# take intersecting genes and combine with neg fdr from DEG
potential_hits = intersect(mus_genes$id, NE_total$X)
# filter hits from DEG df to look at lfc
DE_genes = NE_total |> 
  filter(X %in% potential_hits) |> 
  arrange(log2FoldChange)
# filter hits from the muscka df to look at padj
potential_mus = muscka_100 |> 
  filter(id %in% potential_hits)
# take the genes of that are in the same direction:
common_genes = DE_genes |> 
  filter(log2FoldChange < 0)
mus_genes = mus_genes |> 
  filter(id %in% common_genes$X)

# 5.0 Create csvs for future references ----
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
write.csv(diff_t1, "t1vt3_gene_list.csv") # csv with non-CD8 confounters removed
write.csv(mus_DE, "mus_DEG_combined_hits.csv") # csv for intersecting genese of DEGlo or DEGhi


# 6.0 Calculating z-equivalents of muscka and using the wald statistic ----
z_equiv = function(p, lfc) {
  sign(lfc) * qnorm(1 - p/2)
}
# test to see if padj or pval better indicator for stat
test_z = DE_genes |> 
  transmute(gene = X,
            stat = stat,
            z_equiv = z_equiv(p = pvalue, lfc = log2FoldChange),
            z_equiv_adj = z_equiv(p = padj, lfc = log2FoldChange)) |> 
  arrange(z_equiv)
# conclusion: pvalue is more indicative of Wald Stat than the padj. Note that wald stat
# has less to do with null hypothesis rejection but more so abt no. of SE of lfc is away from 0 (or uncertainty).
z_mus = mus_genes |> 
  transmute(gene = id,
            z_equiv = z_equiv(p = neg.p.value, neg.lfc)) |> 
  arrange(z_equiv)



