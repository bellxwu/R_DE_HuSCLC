## Description: Combining Muscka hits with "infiltration-inducing" DEGs. 
# Goal is to see if combining these hits together and feeding into GSEA will lead to any changes in pathways.
# Used 674 genes from NE-high
# Used 3544 genes from NE-low 
## Author: Bell Wu
## Date created: 2025.09.11

# Changelog:
# 2025.09.18 - Added selection of genes from NE-low also. 

# 1.0 setwd and load dfs -------------------------------------------------------
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
DEGs_hi = read.csv("NE_high_HotvsCold.csv")
dim(DEGs_hi)
head(DEGs_hi)

setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
muscka = read.csv("B6_IgGvsNRG_IgG_0.2.csv")

# 2.0 GSEA of DEGs -----
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)

# create gene list
genelist = DEGs_hi$stat |> 
  setNames(DEGs_hi$X) |> 
  sort(decreasing = TRUE)

# GSEA of NE-high df
GSEA_hi = gseGO(geneList = genelist,
                ont = "BP",
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                pvalueCutoff = 0.2)
view(GSEA_hi)

# 3.0 create Z-equivalent for muscka -----------------------------------------
# function for calculating z-equiv
z_equiv = function(p, lfc) {
  sign(lfc) * qnorm(1 - p/2)
}
# create z-equivalent
z_muscka = muscka |> 
  transmute(gene = id,
            Z_equiv = z_equiv(p = neg.p.value, lfc = neg.lfc))
# combining lists
Z_all = DEGs_hi |> 
  dplyr::select(c("X", "stat"))
colnames(Z_all) = c("gene", "Z_equiv")
Z_all = rbind(Z_all, z_muscka)

# find duplicates:
dup_genes = hi_genelist[duplicated(names(hi_genelist))]
TAP1 = names(hi_genelist) %in% "TAP1"
TAP1 = hi_genelist[TAP1]

# Combined GSEA analysis with GO
hi_genelist = Z_all$Z_equiv |> 
  setNames(Z_all$gene) |> 
  sort(decreasing = TRUE)
head(hi_genelist)

# GSEA analysis
pathways = gseGO(geneList = hi_genelist,
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 pvalueCutoff = 0.2)
head(pathways@result)
view(pathways)

# 4.0 Comparison of previous pathways -----------------------------------------------
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

# 5.0 Visualizing pathways ------------------------------------------------------
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
ggsave("NE-high_only__Delta_NES.png",
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
  ggtitle("GSEA of NE-high Hot vs Cold DEGs") +
  theme_few()
ggsave("NE-high_only_DEG_NES.png",
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
ggsave("NE-high_only_DEGplusMusCKA_NES.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = Mus,
       dpi = 300,
       height = 5,
       width = 10)
# Conclusions: addition of muscka hits to DEG of only NE-low did not result in any significant changes.
# Likely as gene list much larger (674 vs 162).






