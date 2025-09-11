## Description: Combining Muscka hits with "infiltration-inducing" DEGs. 
# Goal is to see if combining these hits together and feeding into GSEA will lead to any changes in pathways.
# Used the 162 genes from both NE-low and NE-high, specifically Wald stat from NE-high.
## Author: Bell Wu
## Date created: 2025.09.05

library(tidyverse)

# 1.0 load dfs and clean data -----
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
muscka = read.csv("B6_IgGvsNRG_IgG_0.2.csv") # FDR < 0.2 muscka

setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
DEG_hi = read.csv("NE_high_sg.csv") # DEGs from NE-high
DEG_lo = read.csv("NE_low_sg.csv") # DEGs from NE-low
DEG_total = rbind(DEG_hi, DEG_lo) # combine for total, note that there will be duplicates

# 2.0 creating z-equivalent as Wald surrogate -----
head(muscka)
# function for calculating z-equiv
z_equiv = function(p, lfc) { 
  sign(lfc) * qnorm(1 - p/2)
}
# create two lists: positive enriched, and negative enriched
# for negative:
z_neg = muscka |> 
  dplyr::transmute(Gene = id,
                   p.value = neg.p.value,
                   padj = neg.fdr,
                   lfc = neg.lfc) |> 
  filter(padj < 0.2)
# for positive: 
z_pos = muscka |> 
  dplyr::transmute(Gene = id,
                   p.value = pos.p.value,
                   padj = pos.fdr,
                   lfc = pos.lfc) |> 
  filter(padj < 0.2)
# combine lists:
Z_muscka = rbind(z_neg, z_pos) |> 
  filter(!grepl("Ctrl", Gene)) # remove ctrls

# calculating df from mageck
Z_muscka = Z_muscka |> 
  dplyr::transmute(Gene = Gene,
                   Z_equiv = z_equiv(p.value, lfc = lfc))
head(Z_muscka)
# testing z-equiv approximation to Wald with DEseq
head(DEG_hi)
Z_DEG_hi = DEG_hi |> 
  dplyr::transmute(Gene = X,
                   Wald = stat,
                   pval = pvalue, 
                   Z_equiv = z_equiv(pvalue, lfc = log2FoldChange))
view(DEG_hi)
head(DEG_hi)
# Not going to floor pvalue, the z_equiv where pval not small similar to Wald

# 3.0 Combining z_equiv to DEG_lo list ---
Z_all = Z_DEG_hi |> 
  dplyr::select(c("Gene", "Wald"))
colnames(Z_all)[2] = "Z_equiv"
dim(Z_all)
head(Z_all)
Z_all = rbind(Z_all, Z_muscka)
dim(Z_all)

# 4.0 Pathway analyses ----
library(clusterProfiler)
library(org.Hs.eg.db)

# create genelist for GO
hi_genelist = Z_all$Z_equiv |> 
  setNames(Z_all$Gene) |> 
  sort(decreasing = TRUE)
head(hi_genelist)

pathways = gseGO(geneList = hi_genelist,
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 pvalueCutoff = 0.5)
head(pathways@result)
view(pathways)

# previous pathways
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
high_res_GO_2 = read.csv("GSEA_GO_NE-high.csv")
view(high_res_GO_2)

# compare NES
DEG_NES = high_res_GO_2 |> 
  dplyr::select(c("ID", "Description", "NES", "p.adjust"))
names(DEG_NES)[-1] <- paste0("DEG_", names(DEG_NES)[-1])
head(DEG_NES)

Mus_NES = pathways@result |> 
  dplyr::select(c("Description", "NES", "p.adjust")) |> 
  rownames_to_column(var = "ID")
names(Mus_NES)[-1] <- paste0("MusCKA_", names(Mus_NES)[-1])
head(Mus_NES)

# join dfs by GO ID
Del_NES = inner_join(DEG_NES, Mus_NES, by = "ID")
head(Del_NES)

# calculate NES deltas and sort
Del_NES = Del_NES |> 
  mutate(Delta_NES = (MusCKA_NES - DEG_NES)) |> 
  dplyr::arrange(desc(Delta_NES))
head(Del_NES)

# 5.0 Visualizing pathways ----
library(ggplot2)
library(ggthemes)

# Visualize DeltaNES
Del_NES = Del_NES |> 
  mutate(significance = ifelse(MusCKA_p.adjust < 0.05, "significant", "non_significant"))
Del = ggplot(Del_NES, mapping = aes(x = Delta_NES, fct_reorder(DEG_Description, Delta_NES), fill = significance)) +
  geom_col() +
  scale_fill_manual(name = "Significance",
                    values = c("significant" = "#E31A1C", "non_significant" = "#1F78B4"),
                    labels = c("significant" = "Significant", "non_significant" = "Non-significant")) +
  ylab(NULL) +
  xlab("Delta NES") +
  ggtitle("Change in NES after inclusion of MusCKA hits") +
  theme_few()
ggsave("DEGhi_MuscKA_Delta_NES.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = Del,
       dpi = 300,
       height = 5,
       width = 10)

# Visualize original NES: DEG_NE_hi
Del_NES = Del_NES |> 
  mutate(significance = ifelse(DEG_p.adjust < 0.05, "significant", "non_significant"))
DEG = ggplot(Del_NES, mapping = aes(x = DEG_NES, fct_reorder(DEG_Description, DEG_NES), fill = significance)) +
  geom_col() +
  scale_fill_manual(name = "Significance",
                    values = c("significant" = "#E31A1C", "non_significant" = "#1F78B4"),
                    labels = c("significant" = "Significant", "non_significant" = "Non-significant")) +
  ylab(NULL) +
  xlab("NES") +
  ggtitle("GSEA of NE-low Hot vs Cold DEGs") +
  theme_few()
ggsave("DEGhi_NES.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = DEG,
       dpi = 300,
       height = 5,
       width = 10)

# Visualize NES after MusCKA
Del_NES = Del_NES |> 
  mutate(significance = ifelse(MusCKA_p.adjust < 0.05, "significant", "non_significant"))
Mus = ggplot(Del_NES, mapping = aes(x = MusCKA_NES, fct_reorder(DEG_Description, MusCKA_NES), fill = significance)) +
  geom_col() +
  scale_fill_manual(name = "Significance",
                    values = c("significant" = "#E31A1C", "non_significant" = "#1F78B4"),
                    labels = c("significant" = "Significant", "non_significant" = "Non-significant")) +
  ylab(NULL) +
  xlab("NES") +
  ggtitle("GSEA of NE-low Hot vs Cold DEGs plus MusCKA hits") +
  theme_few()
ggsave("DEGhi_plus_MusCKA_NES.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = Mus,
       dpi = 300,
       height = 5,
       width = 10)

# 6.0 Identifying leading genes ----
library(enrichplot)
library(ggrepel)
# gseaplot2 of top pathways
# filter for only significant pathways
GO_IDs = Del_NES |> 
  filter(MusCKA_p.adjust < 0.05)
GO_IDs = pathways |> 
  filter(ID %in% GO_IDs$ID)
head(GO_IDs)
dim(GO_IDs)

# Recall the muscka hits
mus_genes = Z_muscka$Gene
head(Z_muscka)
head(GO_IDs$core_enrichment)
core_g = GO_IDs$core_enrichment
# identify leading gene matches from muscka in pathways
results = NULL
for (i in seq_along(core_g)) {
  genes_in_string = strsplit(core_g[i], "/")[[1]]
  matches = intersect(mus_genes, genes_in_string)
  count = length(matches)
  results[[i]] = list(
    description = GO_IDs[ ,2][i],
    # string = core_g[i],
    matches = matches,
    count = count
  )
}
results