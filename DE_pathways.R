## Description: Pathway analyses of DE hits 
## Author: Bell Wu
## Date created: 2025.09.05

BiocManager::install("clusterProfiler")
BiocManager::install("ReactomePA")

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ReactomePA)

# 1.0 setwd and load dfs -----
setwd("~/R_programming/R_DE_HuSCLC/")

# df for DEs in NE-high/ow
NE_high_DE <- read.csv("NE_high_HotvsCold.csv")
NE_low_DE <-read.csv("NE_low_HotvsCold.csv")
dim(NE_high_DE)
dim(NE_low_DE)
head(NE_high_DE)

# 2.0 GSEA-----
# 2.1 test enrichment with gseGO
# order genelist
high_genelist <- NE_high_DE$stat |> 
  setNames(NE_high_DE$X) |> 
  sort(decreasing = TRUE)

# GO pathway
high_res <- gseGO(
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

head(high_res_react)
view(high_res_react)


