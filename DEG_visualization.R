## Description: R script for visualizing genes after differential analysis of selected SCLC cell-lines
## Author: Bell Wu
## Last updated: 2025.02.27

# 1.1: Filtering data ----------------------------------------------------------------
# read csv and clean data
pvals <- read.csv("SCLC_CCLE_Cold-D_vs_Cold-E.csv") 
head(pvals)
row.names(pvals) = pvals$X
pvals <- pvals[, -1]

# to summarize the hits: 
library(dplyr)
pvals0.05 <- pvals %>%
  filter(padj < 0.05) # filter for adj pvals < 0.05
sum(pvals$padj < 0.05, na.rm=TRUE) # find sum total of all pvals < 0.05
dim(pvals0.05) # check dimensions to make sure it matches 

# identify how many genes with significant fold change > a set parameter
upGenes <- pvals0.05 %>%
  filter(log2FoldChange > 1)
dim(upGenes)
downGenes <- pvals0.05 %>%
  filter(log2FoldChange < 1)
dim(downGenes)

# 1.2 Creating a Volcano plot ---------------------------------------------------
# use EnhancedVolcano package to plot
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

EnhVolplot <- EnhancedVolcano(pvals,
                lab = rownames(pvals),
                x = 'log2FoldChange',
                y = 'padj',
                title = "Cold-Excluded vs Cold-Desert",
                selectLab = c('GENE1','GENE2','GENE3',
                              'GENE4','GENE5','GENE6','GENE7','GENE8'),
                subtitleLabSize = 10,
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 1,
                legendPosition = "right",
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1,
                max.overlaps = 10,
                boxedLabels = TRUE
                )
# generate PDF of the plot
pdf(file = "Volcano_plot_DEGs_Cold-E_vs_Cold-D.pdf", width = 10, height = 8)
grid::grid.newpage()
grid::grid.draw(EnhVolplot)
dev.off()

# plot for shrinkLFC, note that there are significant genes but with no 
# fold change so will not continue using this function until further understanding
EnhancedVolcano(DEG_DvE_LFCapel_ordered,
                lab = rownames(DEG_DvE_LFCapel_ordered),
                x = 'log2FoldChange',
                y = 'padj')


