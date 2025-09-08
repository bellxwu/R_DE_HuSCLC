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

# 1.0 setwd and load dfs -----
setwd("~/R_programming/R_DE_HuSCLC/")

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
high_res_GO <- gseGO(
  geneList = high_genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
  )
view(high_res)

# 2.2 ReactomeDB pathway
# convert gene list to EntrezID
eID <- bitr(names(high_genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
high_genelist <- setNames(high_genelist, eID$ENTREZID) # set named vector to entrezID
# calculate 
high_res_react <- gsePathway(
  geneList = high_genelist,
  organism = "human",
  pvalueCutoff = 0.5
)
view(high_res_react)
# Conclusion: in general significant pathways not detected from this, genes that 
# show up from reactomedb did not match any specific pathway curated by reactomeDB
# See companion .md file for in-depth notes.

# 3.0 GSEA of NE-low: -----
# note that NE-low has more genes. NE-low contains 3455 genes.
# 3.0.1 visualize test statistics: 
low_test_stats = data.frame(Gene = NE_low_DE$X,
                            NE_low_stat = NE_low_DE$stat)
high_test_stats = data.frame(Gene = NE_high_DE$X,
                             NE_high = NE_high_DE$stat)
head(test_stats)
# visualize distribution of test statistics for NE-low
low = ggplot(low_test_stats, mapping = aes(x = NE_low_stat)) +
  geom_density(alpha = 0.5)
# visualize for NE-high:
high = ggplot(high_test_stats, mapping = aes(x = NE_high)) +
  geom_density()
low+high

# 3.1 test enrichment with gseGO
# order genelist
low_genelist <- NE_low_DE$stat |> 
  setNames(NE_low_DE$X) |> 
  sort(decreasing = TRUE)

# GO pathway
low_res_GO <- gseGO(
  geneList = low_genelist,
  ont = "BP",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  pvalueCutoff = 0.5
)
view(low_res_GO)

# 3.2 ReactomeDB pathway
# convert gene list to EntrezID
eID <- bitr(names(low_genelist), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
low_genelist <- setNames(low_genelist, eID$ENTREZID) # set named vector to entrezID
# calculate 
low_res_react <- gsePathway(
  geneList = low_genelist,
  organism = "human",
  pvalueCutoff = 0.05
)
view(low_res_react)

# 4.0 Exporting results: write .csv for each -----
# NE-high:
str(high_res)
str(high_res_react)
write.csv(high_res_GO@result, "GSEA_GO_NE-high.csv") # GO
write.csv(high_res_react@result, "GSEA_Reactome_NE-high.csv") # ReactomeDB
# Reactome
write.csv(low_res_GO@result, "GSEA_GO_NE-low.csv") # GO
write.csv(low_res_react@result, "GSEA_Reactome_NE-low.csv") # ReactomeDB



