# Description: Probing TCGA data for normal lung gene expression
# Author: Bell Wu
# Date created: 2025.09.24

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(ggpubr)

# 1.0 load files and read .csv ----------------------------------------
setwd("~/R_programming/R_DE_HuSCLC/Analysis_TCGA/")
TCGA_df = read.table(gzfile("TCGA.LUAD.sampleMap_HiSeqV2.gz"), header = TRUE)
colnames(TCGA_df)[1] = "Gene" 
TCGA_df[1:3, 1:3]

setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
George = read.csv("George_nodup.csv") # data is in FPKM
colnames(George)[1] = "Gene"
head(George)

# load MakeFacetsV2 function
source("~/R_programming/R_analyses_scripts/MakeFacetsV2.R")

# 2.0 Convert FPKM/RPKM to TPM -------------------------------------------------
# 2.1 function to convert; requires genes to be labeled genes ------------------
convert_TPM = function(df) {
  if (!any(colnames(df) == "Gene")) {
    print("Ensure gene column is labeled Gene") # ensure column properly labeled
  } else {
  genes = df |> 
    dplyr::select(Gene) # select out gene names
  df = df |> 
    select(-Gene) # remove gene name from column
  tpm = unlist(apply(df, 2, function(x) {(x / sum(x)) * 1e6})) |> # calculate TPM
    data.frame()
  row.names(tpm) = genes[ , 1]
  }
  return(tpm)
}
# 2.2 Convert datasets to TPM --------------------------------------------------
George_tpm = convert_TPM(George)
TCGA_tpm = convert_TPM(TCGA_df)

dim(George_tpm)
dim(TCGA_tpm)

# double check TPM was successful, should all add up to 1e6
colSums(TCGA_tpm[ , 1:5])
colSums(George_tpm[ , 2:5])

# 3.0 Select only the normal lung samples from TCGA data -----------------------
# Goal is to identify which sampleID relates to normal lung
# TCGA uses identified with "11" to indicate normal adjacent
TCGA_norm = TCGA_tpm[grepl("11$", colnames(TCGA_tpm))] # select sample IDs with 11 at end
# combine George and TCGA
TCGA_norm = rownames_to_column(TCGA_norm, var = "Gene")
George_tpm = rownames_to_column(George_tpm, var = "Gene")
TCGA_George = inner_join(TCGA_norm, George_tpm, by = "Gene")
TCGA_George[1:5,1:5]

# 4.0 Select out the genes of interests (hits undergoing validation): ----------
candidate = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5")
cand_tpm = TCGA_George |> 
  filter(Gene %in% candidate)
# 4.1 Plot out distribution plots ----------------------------------------------
df_long = cand_tpm |> 
  pivot_longer(cols = -Gene,
               names_to = "ID",
               values_to = "TPM")
# Add metadata for sample type
df_long = df_long |> 
  mutate(Sample = case_when(
    str_detect(ID, "TCGA") ~ "Normal"))
df_long$Sample[is.na(df_long$Sample)] = "SCLC"
# add factor levels so smaller sample counts (Normal) will be in front
df_long$Sample = factor(df_long$Sample, levels = ordered(c("SCLC", "Normal")))
# 4.2 plot histogram: ----------------------------------------------------------
# use facets function to plot histogram
hist = MakeFacets(input_df = df_long,
           long = TRUE,
           x = "Sample", # labels for x
           y = "TPM", # labels for x
           fill = "Sample", # variable you want to add colour to
           plot_type = "hist",
           intervals = 11,
           add_stat = FALSE,
           facet_var = "Gene")
# plot histogram: adjust colours:
colours = scale_fill_manual(values = c("Normal" = "#FF0000", "SCLC" = "#5BBCD6"),
                    name = "Sample Type", 
                    labels = c("Normal" = "Normal (n=60)",
                               "SCLC" = "SCLC (n=82)"))
# add labels:
hist = hist +
  colours +
  xlab("TPM Expression bin") +
  ylab("Counts") +
  ggtitle("Histogram of hits in Normal vs SCLC")
# save plot
ggsave("NormalvsSCLC_select_hits_hist.png",
       plot = hist,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 450/90,
       width = 525/90,
       units = "in")
# 4.3 plot boxplots: comparison of medians -------------------------------------
# change order of plot so Normal is first
df_long$Sample = factor(df_long$Sample, levels = ordered(c("Normal", "SCLC")))
# plot box plot:
box = MakeFacets(input_df = df_long,
                 long = TRUE,
                 x = "Sample", # labels for x
                 y = "TPM", # labels for y
                 fill = "Sample", # variable you want to add colour to
                 plot_type = "box",
                 add_stat = TRUE,
                 facet_var = "Gene")
box = box + colours
#save plot
ggsave("NormalvsSCLC_select_hits.png",
       plot = box,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 400/90,
       width = 400/90)

# 5.0 Select out the potential genes (hits from screen <0.2): ------------------
candidate = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6")
cand_tpm = TCGA_George |> 
  filter(Gene %in% candidate)
# 5.1 Plot out distribution plots ----------------------------------------------
df_long = cand_tpm |> 
  pivot_longer(cols = -Gene,
               names_to = "ID",
               values_to = "TPM")
# Add metadata for sample type
df_long = df_long |> 
  mutate(Sample = case_when(
    str_detect(ID, "TCGA") ~ "Normal"))
df_long$Sample[is.na(df_long$Sample)] = "SCLC"
# add factor levels so smaller sample counts (Normal) will be in front
df_long$Sample = factor(df_long$Sample, levels = ordered(c("SCLC", "Normal")))
# 5.2 plot histogram: ----------------------------------------------------------
hist = MakeFacets(input_df = df_long,
                  long = TRUE,
                  x = "Sample", # labels for x
                  y = "TPM", # labels for x
                  fill = "Sample", # variable you want to add colour to
                  plot_type = "hist",
                  intervals = 11,
                  add_stat = FALSE,
                  facet_var = "Gene")
# plot histogram: adjust colours:
colours = scale_fill_manual(values = c("Normal" = "#FF0000", "SCLC" = "#5BBCD6"),
                            name = "Sample Type", 
                            labels = c("Normal" = "Normal (n=60)",
                                       "SCLC" = "SCLC (n=82)"))
# add labels:
hist = hist +
  colours +
  xlab("TPM Expression bin") +
  ylab("Counts") +
  ggtitle("Histogram of hits in Normal vs SCLC")
# save plot
ggsave("NormalvsSCLC_0.2_hits_hist.png",
       plot = hist,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 450/90,
       width = 525/90,
       units = "in")
# 5.3 plot boxplots: comparison of medians -------------------------------------
# change order of plot so Normal is first
df_long$Sample = factor(df_long$Sample, levels = ordered(c("Normal", "SCLC")))
# plot boxplots
box = MakeFacets(input_df = df_long,
                 long = TRUE,
                 x = "Sample", # labels for x
                 y = "TPM", # labels for y
                 fill = "Sample", # variable you want to add colour to
                 plot_type = "box",
                 add_stat = TRUE,
                 facet_var = "Gene")
box = box + colours
#save plot
ggsave("NormalvsSCLC_select_0.2_hits.png",
       plot = box,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 400/90,
       width = 400/90)

# 6.0 Select out genes from NE-low or NE-high ----------------------------------
candidate = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8", "GENE9", "GENE10")
cand_tpm = TCGA_George |> 
  filter(Gene %in% candidate)
# 6.1 Plot out distribution plots ----------------------------------------------
df_long = cand_tpm |> 
  pivot_longer(cols = -Gene,
               names_to = "ID",
               values_to = "TPM")
# Add metadata for sample type
df_long = df_long |> 
  mutate(Sample = case_when(
    str_detect(ID, "TCGA") ~ "Normal"))
df_long$Sample[is.na(df_long$Sample)] = "SCLC"
# add factor levels so smaller sample counts (Normal) will be in front
df_long$Sample = factor(df_long$Sample, levels = ordered(c("SCLC", "Normal")))
# 6.0 plot histogram: ----------------------------------------------------------
# use facets function to plot histogram
hist = MakeFacets(input_df = df_long,
                  long = TRUE,
                  x = "Sample", # labels for x
                  y = "TPM", # labels for x
                  fill = "Sample", # variable you want to add colour to
                  plot_type = "hist",
                  intervals = 11,
                  add_stat = FALSE,
                  facet_var = "Gene")
# plot histogram: adjust colours:
colours = scale_fill_manual(values = c("Normal" = "#FF0000", "SCLC" = "#5BBCD6"),
                            name = "Sample Type", 
                            labels = c("Normal" = "Normal (n=60)",
                                       "SCLC" = "SCLC (n=82)"))
# add labels:
hist = hist +
  colours +
  xlab("TPM Expression bin") +
  ylab("Counts") +
  ggtitle("Histogram of hits in Normal vs SCLC")
# save plot:
ggsave("NormalvsSCLC_DEGplusHits_hist.png",
       plot = hist,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 450/90,
       width = 650/90,
       units = "in")
# 6.3 plot boxplots: comparison of medians -------------------------------------
# change order of plot so Normal is first
df_long$Sample = factor(df_long$Sample, levels = ordered(c("Normal", "SCLC")))
# plot boxplots
box = MakeFacets(input_df = df_long,
                 long = TRUE,
                 x = "Sample", # labels for x
                 y = "TPM", # labels for y
                 fill = "Sample", # variable you want to add colour to
                 plot_type = "box",
                 add_stat = TRUE,
                 facet_var = "Gene")
box = box + colours
#save plot
ggsave("NormalvsSCLC_DEGplusHits_box.png",
       plot = box,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 500/90,
       width = 550/90)

# 7.0 plots for ver2 candidates-------------------------------------------------
# This section is for selection of ver2 candidates to be selected again after more defined criteria filters
candidate = c("MMP11", "DDX23", "TAF7", "SETD1B")
cand_tpm = TCGA_George |> 
  filter(Gene %in% candidate)
# 7.1 Plot out distribution plots ----------------------------------------------
df_long = cand_tpm |> 
  pivot_longer(cols = -Gene,
               names_to = "ID",
               values_to = "TPM")
# Add metadata for sample type
df_long = df_long |> 
  mutate(Sample = case_when(
    str_detect(ID, "TCGA") ~ "Normal"))
df_long$Sample[is.na(df_long$Sample)] = "SCLC"
# add factor levels so smaller sample counts (Normal) will be in front
df_long$Sample = factor(df_long$Sample, levels = ordered(c("SCLC", "Normal")))
# 7.2 plot histogram: ----------------------------------------------------------
# use facets function to plot histogram
hist = MakeFacets(input_df = df_long,
                  long = TRUE,
                  x = "Sample", # labels for x
                  y = "TPM", # labels for x
                  fill = "Sample", # variable you want to add colour to
                  plot_type = "hist",
                  intervals = 11,
                  add_stat = FALSE,
                  facet_var = "Gene")
# plot histogram: adjust colours:
colours = scale_fill_manual(values = c("Normal" = "#FF0000", "SCLC" = "#5BBCD6"),
                            name = "Sample Type", 
                            labels = c("Normal" = "Normal (n=60)",
                                       "SCLC" = "SCLC (n=82)"))
# add labels:
hist = hist +
  colours +
  xlab("TPM Expression bin") +
  ylab("Counts") +
  ggtitle("Histogram of hits in Normal vs SCLC")
# save plot:
ggsave("NormalvsSCLC_DEGplusHits_hist_candidates-ver2.png",
       plot = hist,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 450/90,
       width = 650/90,
       units = "in")
# 7.3 plot boxplots: comparison of medians -------------------------------------
# change order of plot so Normal is first
df_long$Sample = factor(df_long$Sample, levels = ordered(c("Normal", "SCLC")))
# plot boxplots
box = MakeFacets(input_df = df_long,
                 long = TRUE,
                 x = "Sample", # labels for x
                 y = "TPM", # labels for y
                 fill = "Sample", # variable you want to add colour to
                 plot_type = "box",
                 add_stat = TRUE,
                 facet_var = "Gene")
box = box + colours
#save plot
ggsave("NormalvsSCLC_DEGplusHits_box_candidates-ver2.png",
       plot = box,
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       dpi = 300,
       height = 500/90,
       width = 550/90)
