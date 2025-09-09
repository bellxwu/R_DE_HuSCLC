## Description: Pathway analyses of DE hits for NE-high and NE-low. Goal of this analysis is to determine if
# pathways shown make sense and can be corroborated with MusCKA hits. 
## Author: Bell Wu
## Date created: 2025.09.05

BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)
library(ggplot2)
library(enrichplot)

# 1.0 setwd and load dfs -----
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")

# df for DEs in NE-high/ow
NE_high_DE <- read.csv("NE_high_HotvsCold.csv")
NE_low_DE <-read.csv("NE_low_HotvsCold.csv")
dim(NE_high_DE)
dim(NE_low_DE)
head(NE_high_DE)

# 2.0 GSEA of NE-high -----
# NE-high gene list contains 674 genes
# 2.1 test enrichment with gseGO
# order genelist
high_genelist <- NE_high_DE$stat |> 
  setNames(NE_high_DE$X) |> 
  sort(decreasing = TRUE)

# GO pathway
high_res_GO_2 <- gseGO(
  geneList = high_genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
  )
view(high_res_GO_2)

# 2.2 ReactomeDB pathway
# convert gene list to EntrezID
eID <- bitr(names(high_genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
high_genelist <- setNames(high_genelist, eID$ENTREZID) # set named vector to entrezID
# calculate 
high_res_react_2 <- gsePathway(
  geneList = high_genelist,
  organism = "human",
  pvalueCutoff = 0.5
)
view(high_res_react_3)
# Conclusion: in general significant pathways not detected from this, genes that 
# show up from reactomedb did not match any specific pathway curated by reactomeDB
# See companion .md file for in-depth notes.

# 3.0 GSEA of NE-low: -----
# note that NE-low has more genes. NE-low contains 3455 genes.
# 3.0.1 visualize test statistics: ----
low_test_stats = data.frame(Gene = NE_low_DE$X,
                            NE_low_stat = NE_low_DE$stat)
high_test_stats = data.frame(Gene = NE_high_DE$X,
                             NE_high = NE_high_DE$stat)
# visualize distribution of test statistics for NE-low
low = ggplot(low_test_stats, mapping = aes(x = NE_low_stat)) +
  geom_density(alpha = 0.5)
# visualize for NE-high:
high = ggplot(high_test_stats, mapping = aes(x = NE_high)) +
  geom_density()
low+high

# 3.1 test enrichment with gseGO ------
# order genelist
low_genelist <- NE_low_DE$stat |> 
  setNames(NE_low_DE$X) |> 
  sort(decreasing = TRUE)

# GO pathway
low_res_GO_3 <- gseGO(
  geneList = low_genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
)
view(low_res_GO_3)

# 3.2 ReactomeDB pathway
# convert gene list to EntrezID
eID <- bitr(names(low_genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
low_genelist <- setNames(low_genelist, eID$ENTREZID) # set named vector to entrezID
# calculate 
low_res_react_3 <- gsePathway(
  geneList = low_genelist,
  organism = "human",
  pvalueCutoff = 0.05
)
view(low_res_react)

# 4.0 Exporting results: write .csv for each -----
# NE-high:
str(high_res_GO_2)
str(high_res_react_2)
write.csv(high_res_GO_2@result, "~/R_programming/R_DE_HuSCLC/Analysis_csvs/GSEA_GO_NE-high.csv") # GO
write.csv(high_res_react_2@result, "~/R_programming/R_DE_HuSCLC/Analysis_csvs/GSEA_Reactome_NE-high.csv") # ReactomeDB
# Reactome
write.csv(low_res_GO_3@result, "~/R_programming/R_DE_HuSCLC/Analysis_csvs/GSEA_GO_NE-low.csv") # GO
write.csv(low_res_react_3@result, "~/R_programming/R_DE_HuSCLC/Analysis_csvs/GSEA_Reactome_NE-low.csv") # ReactomeDB

# 5.0 Combining DEG list from NE-high and NE-low: -----
# Use intersecting genes and need to make two tests of pathways based on different Wald

# 5.1 gene-list chosen from NE-high:
# gene list for NE-high
high_genelist <- NE_high_DE$stat |> 
  setNames(NE_high_DE$X) |> 
  sort(decreasing = TRUE)
# gene list for NE-low
low_genelist <- NE_low_DE$stat |> 
  setNames(NE_low_DE$X) |> 
  sort(decreasing = TRUE)
# find intersecting genes
same_genes = intersect(names(high_genelist), names(low_genelist)) # intersecting genes
same_genes_high = high_genelist[names(high_genelist) %in% same_genes] 
# pathway analysis
high_res_GO_5 <- gseGO(
  geneList = same_genes_high,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
)
view(high_res_GO_5)
# shows significant pathway

# 5.2 gene-list on NE-low
same_genes_low = low_genelist[names(low_genelist) %in% same_genes]
# pathway analysis
low_res_GO_5 <- gseGO(
  geneList = same_genes_low,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
)
view(low_res_GO_5)
# no enriched pathways at all...
# conclusions: Pathways shown dependent on which DEG Wald-statistic was taken from.
# Directionality is different for intersecting genes thus pathways shown will be different.

# 5.3 Test pathways with reactomeDB
# using NE-high Walds
eID <- bitr(names(same_genes_high), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
same_genes_high <- setNames(same_genes_high, eID$ENTREZID) # set named vector to entrezID
high_res_react <- gsePathway(
  geneList = same_genes_high,
  organism = "human",
  pvalueCutoff = 0.2
)
view(high_res_react)

# using NE-low Walds
# convert gene list to EntrezID
eID <- bitr(names(same_genes_low), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
same_genes_low <- setNames(same_genes_low, eID$ENTREZID) # set named vector to entrezID
low_res_react <- gsePathway(
  geneList = same_genes_low,
  organism = "human",
  pvalueCutoff = 0.05
)
view(low_res_react)
# few pathways shown, also pathways shown more general

# 6.0 Visualizing GSEA pathways: ----
# bar plot of enrichment scores
# from results of GO in 2.0
plot_2 <- ggplot(high_res_GO_2, showCategory=10, aes(NES, fct_reorder(Description, NES), fill = qvalue)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#E31A1C", "#1F78B4" ,"#1F78B4"),
                       guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Hot vs Cold DEG pathway")
# save plot
ggsave("NE-high_DEG_pathway.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       plot = plot_2)

# from results of GO in 5.0
plot_5 <- ggplot(high_res_GO_2, showCategory=10, aes(NES, fct_reorder(Description, NES), fill = qvalue)) +
  geom_col() +
  scale_fill_gradientn(colours=c("#E31A1C", "#1F78B4" ,"#1F78B4"),
                       guide=guide_colorbar(reverse=TRUE)) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("Hot vs Cold DEG pathway")
# save plot
ggsave("NE-high_combined_DEG_pathway.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       plot = plot_5)





