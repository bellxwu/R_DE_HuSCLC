## Description: Compiling the datasets from CCLE and editing to include the 
#cell-lines of interest
## Author: Bell Wu
## Last updated: 2025.02.04

# Set working directory:
getwd()
setwd("/Users/bellwu/R_programming/R_DE_HuSCLC")
setwd("~/R_programming/R_DE_HuSCLC/")

# read csv files
SCLC_CCLE <- read.csv("SCLC_CCLE.csv")
subset_SCLC <- head(SCLC_CCLE) # observe layout

library(tidyr) 
#code with smaller subset
subset_SCLC <- separate(subset_SCLC, X, into = c("Gene", "GeneID"), 
                        sep = "\\.\\.")
row.names(subset_SCLC) = subset_SCLC$Gene
subset_SCLC <- subset_SCLC[ , c(-1,-2)]

# apply code to larger subset
SCLC_CCLE <- separate(SCLC_CCLE, X, into = c("Gene", "GeneID"), 
                        sep = "\\.\\.")
# NAs are detected, will need to identify them
incomplete_rows <- SCLC_CCLE[!complete.cases(SCLC_CCLE), ] # identify which row has the NA
dim(SCLC_CCLE) # check dimensions
SCLC_CCLE <- SCLC_CCLE[complete.cases(SCLC_CCLE), ]
dim(SCLC_CCLE) # check dimensions again to see if it worked
# Great success! - Borat

row.names(SCLC_CCLE) = SCLC_CCLE$Gene
SCLC_CCLE <- SCLC_CCLE[ , c(-1,-2)]

library(dplyr)
selected_lines <- SCLC_CCLE %>%
  select(contains("H841"), contains("SHP77"), contains("H82"), contains("H446"), 
         contains("SBC5"), contains("H1694")
         )
# after running this script, will have a list of cell-lines of interest



