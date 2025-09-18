# Description: Validation by prevalence of nominated sgRNA hits on George dataset.
# Goal is to see if the gene hits from muscka and DEGs are high in prevalence.
# Author: Bell Wu
# Date created: 2025.09.17

library(tidyverse)
library(ggthemes)
library(RColorBrewer)

# 1.0 set wd and load all dfs ----
# Load George dataset 
setwd("~/R_programming/R_RPP_macrophage_adenosine/")
George = read.csv("George_nodup.csv") # data is in FPKM
colnames(George)[1] = "Gene"
head(George)

# Load muscka hits
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
# load < 0.2 FDR hits
muscka_0.2 = read.csv("B6_IgGvsNRG_IgG_0.2.csv") # FDR < 0.2 muscka
muscka_0.2 = muscka_0.2 |> 
  filter(!grepl("Ctrl", id)) # remove Ctrl ID's
# load top 100 neg hits
muscka_100 = read.csv("B6_IgGvsNRG_IgG_neg100.csv")
muscka_100 = muscka_100 |> 
  filter(!grepl("Ctrl", id)) # remove Ctrl ID's
dim(muscka_100)

# Load DEGs
setwd("~/R_programming/R_DE_HuSCLC/Analysis_csvs/")
DEG_lo = read.csv("NE_low_sg.csv") # same gene DEGs w NE-low Wald Stat

# Load CCLE data
setwd("~/R_programming/R_CCLE/")
CCLE = read.csv("CLEAN_SCLC_CCLE.csv")
colnames(CCLE)[1] = "Gene"

# 2.0 select candidate hits from dataset ----
# 2.1 for candidate hits
candidates = c("PEO1", "DDX23", "TAF7", "NRK", "PECAM1", "SLC6A11", "FBXO11", "KDM1A")
# from George data set
George_hits = George |> 
  filter(Gene %in% candidates)
dim(George_hits)
# from CCLE
CCLE_hits = CCLE |> 
  filter(Gene %in% candidates)
dim(CCLE_hits)

# alternate gene names
peo1 = c("TWINKLE", "IOSCA", "MTDPS7", "PRLTS5", "ATXN8", "PEOA3", "SANDO", "SCA8", "PEO1", "TWINL", "PEO", "FLJ21832")
pecam1 = c("PECA1", "CD31", "PECAM1", "PECAM", "ENDOCAM")

CCLE_hits = CCLE |> 
  filter(Gene %in% pecam1)

# 2.2 for top 100 neg hits
candidates_100 = muscka_100$id
CCLE_hits = CCLE |> 
  filter(Gene %in% candidates_100)

# 3.0 Average expression levels -----
# for chosen hits
CCLE_means = CCLE_hits |> 
  transmute(Gene = Gene,
            Row_means = rowMeans(CCLE_hits[-1]))

# for all neg hits
CCLE_means = CCLE_hits |> 
  transmute(Gene = Gene,
            Row_means = rowMeans(CCLE_hits[-1]))

# 4.0 plot candidate hit levels across patients -----
# Purpose of this section is to plot the candidate hits across the different datasets

# 4.1 Plot out the candidate hits for CCLE -------
s_cols = colnames(CCLE_hits)[-1] # take samples 
d_candidate = CCLE_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# plot density for the candidate genes selected previously in CCLE
d_p = ggplot(data = d_candidate, mapping = aes(x = expression, fill = Gene)) +
  geom_density(alpha = 0.75,
               linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  # scale_color_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Density plots of RNA \nexpression in CCLE")
# save plot
ggsave("density_CCLE_RNA_hits.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = d_p,
       dpi = 300,
       height = 3,
       width = 4)

# 4.2 plot from George_2015 dataset -------
s_cols = colnames(George_hits)[-1] # take samples 
d_candidate = George_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# set bins
breaks = c(seq(0, 70, by = 10), Inf)
d_candidate = d_candidate |> 
  mutate(bin_expression = cut(expression, breaks = breaks, right = FALSE,
                              labels = c("0–10", "10–20", "20–30", "30–40", 
                                         "40–50", "50–60", "60–70", "70+")))
# plot histogram for the candidate genes selected previously in George
h_p = ggplot(data = d_candidate, mapping = aes(x = bin_expression, fill = Gene)) +
  geom_bar(position = "dodge",
           color = "black",
           linewidth = 0.25,
           stat = "count") + 
  scale_fill_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Histogram of RNA expression \nin George 2015") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# save plot
ggsave("histogram_George_RNA_hits.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = h_p,
       dpi = 300,
       height = 3,
       width = 4.5)

# 5.0 Analysis of potential hits (DEG + muscka) ----
# Purpose is to find intersecting genes from DEG-hi or DEG-lo that are similar to muscka hits
setwd("/Users/bellwu/R_programming/R_DE_HuSCLC/Analysis_csvs")
potential_hits = read.csv("mus_DEG_combined_hits.csv") # read csv for combined hits
head(potential_hits)
candidates = potential_hits$id # take gene names

# 5.1 Pull candidate genes from dataset ------
# from George dataset
George_hits = George |> 
  filter(Gene %in% candidates, !Gene %in% c("NRK", "TAP1", "SLC6A11"))
dim(George_hits)
# from CCLE
CCLE_hits = CCLE |> 
  filter(Gene %in% candidates, !Gene %in% c("NRK", "TAP1", "SLC6A11"))
dim(CCLE_hits)

# 5.2 plot density/histogram plots from dataset ----
s_cols = colnames(CCLE_hits)[-1] # take samples 
d_candidate = CCLE_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# plot density for the candidate genes selected previously in CCLE
d_p = ggplot(data = d_candidate, mapping = aes(x = expression, fill = Gene)) +
  geom_density(alpha = 0.75,
               linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  # scale_color_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Density plots of RNA \nexpression in CCLE")
# save plot
ggsave("density_CCLE_RNA_hits_NE_total.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = d_p,
       dpi = 300,
       height = 3,
       width = 4)

# 5.2 plot from George_2015 dataset -------
s_cols = colnames(George_hits)[-1] # take samples 
d_candidate = George_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# set bins
breaks = c(seq(0, 70, by = 10), Inf)
d_candidate = d_candidate |> 
  mutate(bin_expression = cut(expression, breaks = breaks, right = FALSE,
                              labels = c("0–10", "10–20", "20–30", "30–40", 
                                         "40–50", "50–60", "60–70", "70+")))
# plot histogram for the candidate genes selected previously in George
h_p = ggplot(data = d_candidate, mapping = aes(x = bin_expression, fill = Gene)) +
  geom_bar(position = "dodge",
           color = "black",
           linewidth = 0.25,
           stat = "count") + 
  scale_fill_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Histogram of RNA expression \nin George 2015") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("histogram_George_RNA_hits_NE_total.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = h_p,
       dpi = 300,
       height = 3,
       width = 4.5)

# 6.0 Analysis of potential hits (MusCKA < 0.2) ------
# Purpose is to pull hits of < 0.2 FDR
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver2_analyses/comparative_analysis/")
candidates = muscka_0.2$id

# 6.1 Pull data from datasets and remove duplicates: ------
# from George dataset
George_hits = George |> 
  filter(Gene %in% candidates, !Gene %in% c("TAP1", "B2M", "HNRNPK", "DDX23", "TICRR", "PECAM1", "AURKB", "TAF7"))
dim(George_hits)
# from CCLE
CCLE_hits = CCLE |> 
  filter(Gene %in% candidates, !Gene %in% c("TAP1", "B2M", "HNRNPK", "DDX23", "TICRR", "PECAM1", "AURKB", "TAF7"))
dim(CCLE_hits)

# 6.2 Plot density from CCLE -------
s_cols = colnames(CCLE_hits)[-1] # take samples 
d_candidate = CCLE_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# plot density for the candidate genes selected previously in CCLE
d_p = ggplot(data = d_candidate, mapping = aes(x = expression, fill = Gene)) +
  geom_density(alpha = 0.75,
               linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  # scale_color_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Density plots of RNA \nexpression in CCLE")
# save plot
ggsave("density_CCLE_RNA_hits_mus0.2.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = d_p,
       dpi = 300,
       height = 3,
       width = 4)

# 6.2 plot from George_2015 dataset -------
s_cols = colnames(George_hits)[-1] # take samples 
d_candidate = George_hits |>  # create a df_long
  pivot_longer(cols = s_cols,
               names_to = "samples",
               values_to = "expression")
# set bins
breaks = c(seq(0, 70, by = 10), Inf)
d_candidate = d_candidate |> 
  mutate(bin_expression = cut(expression, breaks = breaks, right = FALSE,
                              labels = c("0–10", "10–20", "20–30", "30–40", 
                                         "40–50", "50–60", "60–70", "70+")))
# plot histogram for the candidate genes selected previously in George
h_p = ggplot(data = d_candidate, mapping = aes(x = bin_expression, fill = Gene)) +
  geom_bar(position = "dodge",
           color = "black",
           linewidth = 0.25,
           stat = "count") + 
  scale_fill_brewer(palette = "Set3") +
  theme_few() +
  ggtitle("Histogram of RNA expression \nin George 2015") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("histogram_George_RNA_hits_mus0.2.png",
       path = "~/R_programming/R_DE_HuSCLC/Analysis_PNGs/",
       plot = h_p,
       dpi = 300,
       height = 3,
       width = 4.5)

