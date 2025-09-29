# Description: QC of B6-IgG vs NRG-IgG. Includes gini indices, % counts, etc.
# Ensure quality of sgRNA before downstream analysis
# Date created: 2025.09.02
# Author: Bell Wu

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggthemes)
library(wesanderson)

# 1.0 prepare dfs and workspace -----
setwd("~/R_programming/R_MuSCK_Library/MusCKA-run3/Ver1_analyses/")

# load data files:
counts_summary <- read.table("MuSCK3_all.countsummary.txt", header = TRUE)
counts <- read.table("MuSCK3_all.count.txt", header = TRUE)
counts_norm <- read.table("MuSCK3_all.count_normalized.txt", header = TRUE)

head(counts_summary)
head(counts)
head(counts_norm)

# remove ctrl sgRNAs
no_ctrl = counts |> 
  filter(!grepl("Ctrl", sgRNA))
dim(counts)
dim(no_ctrl)

# 2.0 create barplot for the average of sgRNAs mapped ----
counts_avg <- counts_summary |> 
  group_by(Label) |> 
  summarise(avg_zerocounts = mean(Zerocounts)) |> 
  mutate(sgrnas_detected = 24622 - avg_zerocounts)

ggplot(counts_avg, mapping = aes(x = reorder(Label, -sgrnas_detected), y = sgrnas_detected, fill = Label)) +
  geom_col() +
  labs(title = "Average sgRNAs detected per sample",
       x = "Sample", y = "sgRNAs detected") +
  scale_fill_manual(values = wes_palette("BottleRocket1")[1:6]) +
  geom_text(aes(label = scales::comma(sgrnas_detected)),   # put the n values
            vjust = -0.3, size = 3.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # rotate x axis

ggsave(
  file = "Avg_sgRNAs.png",
  plot = last_plot(),
  path = "~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/sgRNA_QC",
  width = 5,
  height = 5,
  dpi = 300
)

# 2.1 create bar plot for percent mapped and average gini indices -----
# 2.1.1 percent reads mapped:-----
map_avg <- counts_summary |> 
  group_by(Label) |> 
  summarise(percent_mapped = mean(Percentage), average_Gini = mean(GiniIndex))

ggplot(map_avg, mapping = aes(x = reorder(Label, -percent_mapped), y = percent_mapped, fill = Label)) +
  geom_col() +
  labs(title = "Percent reads mapped sgRNAs per group",
       x = "Sample", y = "sgRNAs Mapped/sgRNAs reads") +
  scale_fill_manual(values = wes_palette("BottleRocket1")[1:6]) +
  geom_text(aes(label = scales::comma(percent_mapped)),   # put the n values
            vjust = -0.3, size = 3.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # rotate x axis
# save plot
ggsave(
  file = "Percent_reads_mapped.png",
  plot = last_plot(),
  path = "~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/sgRNA_QC",
  width = 5,
  height = 5,
  dpi = 300
)

# 2.1.2 Gini indicies: -----
ggplot(map_avg, mapping = aes(x = reorder(Label, -average_Gini), y = average_Gini, fill = Label)) +
  geom_col() +
  labs(title = "Gini indices per sample",
       x = "Sample", y = "Gini index") +
  scale_fill_manual(values = wes_palette("BottleRocket1")[1:6]) +
  geom_text(aes(label = scales::comma(average_Gini)),   # put the n values
            vjust = -0.3, size = 3.5) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) # rotate x axis
# save plot
ggsave(
  file = "Gini_indices.png",
  plot = last_plot(),
  path = "~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/sgRNA_QC",
  width = 5,
  height = 5,
  dpi = 300
)

# 3.0 Density plots of genes ----
head(counts_norm)
# create a long df
norm_long <- counts_norm |> 
  dplyr::select(-Gene) |> 
  pivot_longer(-c("sgRNA"),
               names_to = "Sample",
               values_to = "Counts")
head(norm_long)
# Create density plot
ggplot(norm_long, aes(x = log10(Counts + 1), color = Sample, fill = Sample)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Density of normalized sgRNA counts",
    x = "log10(normalized counts + 1)",
    y = "Density"
  ) +
  theme_bw()
# save plot
ggsave(
  file = "Sample_densities.png",
  plot = last_plot(),
  path = "~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/sgRNA_QC",
  width = 5,
  height = 5,
  dpi = 300
)

# remove the aCD8_T1 data from future analysis
norm_long <- counts_norm |> 
  dplyr::select(-c("Gene", "B6_aCD8_T1")) |> 
  pivot_longer(-c("sgRNA"),
               names_to = "Sample",
               values_to = "Counts")
# Create density plot
ggplot(norm_long, aes(x = log10(Counts + 1), color = Sample, fill = Sample)) +
  geom_density(alpha = 0.3) +
  labs(
    title = "Density of normalized sgRNA counts",
    x = "log10(normalized counts + 1)",
    y = "Density"
  ) +
  theme_bw()
# save plot
ggsave(
  file = "Sample_densities_no_aCD8_T1.png",
  plot = last_plot(),
  path = "~/Desktop/Lok Lab/*HuMice Project - CRISPR screen/MusCKA-RUN3/sgRNA_QC",
  width = 5,
  height = 5,
  dpi = 300)

