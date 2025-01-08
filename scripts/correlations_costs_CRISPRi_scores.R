
setwd("~/Documents/CRISPRi")

library(ggplot2)
library(xlsx)
library(dplyr)
library(tidyr)
library(paletteer)
library(stringr)
library(plotly)
library(htmlwidgets)
library(tidyverse)
library(writexl)
library(ggpubr)
library(car)

# Import all data (2nd and 1st demultiplexing)

median_log2FC <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene_definitive.xlsx", 1)
metadata <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene_definitive.xlsx", 2)

annot_info <- read.xlsx("/home/jorge/Documents/CRISPRi/pOXA48_annot.xlsx", 1)
new_annot_info <- read.xlsx("/home/jorge/Documents/CRISPRi/new_annot_pOXA48.xlsx", 1)

medianlog_df <- median_log2FC %>%
  pivot_longer(!gene, names_to = "sample_ID", values_to = "median_log2FC")

sample_metadata <- metadata %>% select(sample_ID, strain, timepoint, treatment)
sample_metadata <- unique( sample_metadata[,c("sample_ID", "strain", "timepoint", "treatment")])

mlog2FC_w_metadata <- merge(medianlog_df, sample_metadata, by = "sample_ID") # Include sample information (strain column)
mlog2FC <- merge(mlog2FC_w_metadata, annot_info, by = "gene") # Include gene annotation information
mlog2FC <- merge(mlog2FC, new_annot_info, by = "gene") # include intergenic regions information

# Include species information in mlog2FC table

strain <- sort(unique(mlog2FC[,"strain"]))
species <- c(rep("E. coli", 2), rep("K. pneumoniae", 2), rep("E. coli", 4), rep("K. pneumoniae", 6))
species_info <- as.data.frame(cbind(strain, species))

mlog2FC <- merge(mlog2FC, species_info, by = "strain")

# Change strain names in metadata to remove the e at the beginning which indicates the treatment already shown in the treatment column

for(i in 1:nrow(`mlog2FC`)) {       # for-loop over rows
  if (str_sub(mlog2FC$strain[i], start = 1, end = length(mlog2FC$strain[i])) == "e") {
    mlog2FC$strain[i] <- str_sub(mlog2FC$strain[i], start = 2)
  }
}

cost_info <- read.xlsx("/home/jorge/Documents/CRISPRi/correlations_costs/costes_curvas_cepas.xlsx", sheetIndex = 2)

mlog2FC <- merge(mlog2FC, cost_info, by = "strain")
mlog2FC$coste_OXA <- as.numeric(mlog2FC$coste_OXA)
mlog2FC$median_log2FC <- as.numeric(mlog2FC$median_log2FC)

# Plot correlation between costs and blaOXA-48 scores

correlations_costs_blaoxa <- mlog2FC %>%
  filter(treatment == "LB+DAPG", timepoint == "t3") %>%
  ggplot(aes(x = coste_OXA, y = median_log2FC, color = species)) +
  geom_point()+
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.1)+
  stat_cor(method = "pearson", size = 1.5)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  #geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 5)+
  facet_wrap(~annot_gene)+
  xlab("Relative cost pOXA-48")+
  ylab("Gene score (Log2 FC)")+
  #ggtitle("DNDJGHEP_00002/blaOXA-48_A score vs pOXA-48 & pFR cost correlation merged") +
  #ylim(-10.5,10.5)+
  #scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(correlations_costs_blaoxa)

# Plot correlation vs blaOXA-48 score (Ec Kpn together)

mlog2FC %>%
  filter(treatment == "LB+DAPG", timepoint == "t2", gene == "blaOXA-48") %>%
  ggplot(aes(x = coste_OXA, y = median_log2FC), col = species) +
  geom_point(aes(color = species))+
  geom_smooth(aes(col = species), method = "lm", se = TRUE, col = "black", alpha = 0.1)+
  stat_cor(aes(col = species),method = "pearson", size = 5)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  #geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 18)+
  facet_wrap(~annot_gene)+
  xlab("Relative cost pOXA-48")+
  ylab("Gene score (Log2 FC)")+
  #ggtitle("DNDJGHEP_00002/blaOXA-48_A score vs pOXA-48 & pFR cost correlation merged") +
  #ylim(-10.5,10.5)+
  #scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

### Check normality

data_stats <- mlog2FC %>%
  filter(treatment == "LB+DAPG", timepoint == "t3", gene == "blaOXA-48")

model <- lm(median_log2FC ~ coste_OXA + species, data = data_stats)

ggqqplot(residuals(model))

shapiro.test(residuals(model)) # Normal data -> Pearson correlation


mlog2FC %>%
  filter(treatment == "LB+DAPG", timepoint == "t3", gene == "blaOXA-48") %>%
  ggplot(aes(x = coste_OXA, y = median_log2FC, color = species)) +
  geom_point(aes(color = species), size = 2.5)+
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1)+
  stat_cor(method = "pearson", size = 5)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  #geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 18)+
  facet_wrap(~annot_gene)+
  xlab("Relative cost pOXA-48")+
  ylab("Gene score (Log2 FC)")+
  #ggtitle("DNDJGHEP_00002/blaOXA-48_A score vs pOXA-48 & pFR cost correlation merged") +
  #ylim(-10.5,10.5)+
  #scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

