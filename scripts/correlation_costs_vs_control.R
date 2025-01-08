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
library(reshape2)

median_log2FC <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_guide_definitive.xlsx", 1)
metadata <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_guide_definitive.xlsx", 2)
annot_info <- read.xlsx("/home/jorge/Documents/CRISPRi/pOXA48_annot.xlsx", 1)
new_annot_info <- read.xlsx("/home/jorge/Documents/CRISPRi/new_annot_pOXA48.xlsx", 1)

cost_info <- read.xlsx("/home/jorge/Documents/CRISPRi/correlations_costs/costes_curvas_cepas.xlsx", 1)

colnames(median_log2FC) <- as.character(colnames(median_log2FC))

medianlog_df <- median_log2FC %>%
  pivot_longer(!c(gene, name, guide), names_to = "sample_ID", values_to = "median_log2FC")

medianlog_df <- as.data.frame(medianlog_df)

sample_metadata <- metadata %>% select(sample_ID, strain, timepoint, treatment)
sample_metadata <- unique( sample_metadata[,c("sample_ID", "strain", "timepoint", "treatment")])

mlog2FC_w_metadata <- merge(medianlog_df, sample_metadata, by = "sample_ID") # Include sample information (strain column)

# Plot distributions of log2FC; but first, filter table and split between control and interference groups

mlog2FC_w_metadata <- mlog2FC_w_metadata %>%
  filter(timepoint == "t3")

mlog2FC_grouped <- transform(mlog2FC_w_metadata, group = ifelse(gene=="control", "control", "interference"))

mlog2FC_grouped <- mlog2FC_grouped %>% 
  filter(treatment == "LB+DAPG")

mlog2FC_grouped <- mlog2FC_grouped %>%
  drop_na()

mlog2FC_grouped %>%
  ggplot(aes(x = median_log2FC, col = group)) +
  geom_density() +
  facet_wrap(~strain) +
  geom_vline(aes(xintercept = mean(median_log2FC), color = group), linetype = "dashed") +
  theme_bw()+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

meantables <- aggregate(x = mlog2FC_grouped$median_log2FC, by = list(mlog2FC_grouped$group, mlog2FC_grouped$strain), FUN = mean)

colnames(meantables) <- c("group", "strain", "mean")

meanscontrol <- meantables %>%
  filter(group == "control")

meansinterference <- meantables %>%
  filter(group == "interference")

meanscontrol$delta <- meansinterference$mean - meanscontrol$mean

corrtable <- merge(meanscontrol, cost_info, by = "strain")

corrtable$species <- c("E. coli", "E. coli", "K. pneumoniae", "K. pneumoniae",
                       rep("E. coli", 4), rep("K. pneumoniae", 6))

model <- lm(coste_OXA ~ delta, data = corrtable)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

corrtable %>%
  ggplot(aes(x = coste_OXA, y = delta), col = species) +
  geom_point(aes(col = species), size = 3) +
  stat_cor(aes(col = species), method = "spearman", label.x.npc = 0.8) +
  geom_smooth(method = "lm", alpha = 0.1, col = "black") +
  ylab("Δlog2FC (CRISPRi - control)") +
  xlab("Relative fitness") +
  theme_bw(base_size = 14)+
  theme(panel.background = element_blank(), panel.grid = element_blank(), 
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

corrtable_coli <- corrtable %>%
  filter(species == "E. coli")

model <- lm(coste_OXA ~ delta, data = corrtable_coli)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

corrtable_kpn <- corrtable %>%
  filter(species == "K. pneumoniae")

model <- lm(coste_OXA ~ delta, data = corrtable_kpn)
ggqqplot(residuals(model))

shapiro.test(residuals(model))

corrtable %>%
  ggplot(aes(x = coste_OXA, y = delta)) +
  geom_point(aes(col = species),size = 3) +
  stat_cor(method = "pearson") +
  #facet_wrap(~species, scales = "free_x") +
  geom_smooth(method = "lm", alpha = 0.1, col = "black") +
  ylab("Δlog2FC (CRISPRi - control)") +
  xlab("Cost pOXA-48") +
  theme_bw(base_size = 14)+
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        strip.text = element_text(face = "italic"),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))
