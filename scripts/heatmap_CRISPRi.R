
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

# Include the new data for the second demultiplexing
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

# Get only timepoint 3

mlog2FC <- mlog2FC %>%
  filter(timepoint == "t3")

mlog2FC %>%
  pivot_wider(names_from = annot_gene, values_from = median_log2FC)

mlog2FC %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_noab <- mlog2FC %>%
  filter(treatment == "LB+DAPG") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_ab <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

mlog2FC %>%
  ggplot(aes(x = factor(annot_gene, levels = c('repA', 'repA/trbC_A', 'repA/trbC_B',
                                               'trbC', 'trbB', 'trbA', 'trbN', 'DNDJGHEP_00001',
                                               'tir (1/2)', 'IS4-family_1', 'IS1-family_1', 'DNDJGHEP_00002',
                                               'DNDJGHEP_00002/blaOXA-48_A', 'DNDJGHEP_00002/blaOXA-48_B',
                                               'blaOXA-48', 'blaOXA-48/lysR_A', 'blaOXA-48/lysR_B', 'lysR',
                                               'lysR/IS4-family_A', 'lysR/IS4-family_B', 'IS4-family_2', 'tir (2/2)',
                                               'pemI', 'pemK', 'ltrA', 'ltrA/umuD_A', 'ltrA/umuD_B', 'umuD', 'umuC', 
                                               'DNDJGHEP_00003', 'DNDJGHEP_00004', 'DNDJGHEP_00005',
                                               'helix-turn-helix transcriptional regulator-1', 'relB',
                                               'relB/DNDJGHEP_00006_A', 'relB/DNDJGHEP_00006_B', 'DNDJGHEP_00006',
                                               'DNDJGHEP_00007', 'DNDJGHEP_00008', 'DNDJGHEP_00008/DNDJGHEP_00009_A',
                                               'DNDJGHEP_00008/DNDJGHEP_00009_B', 'DNDJGHEP_00009', 'DNDJGHEP_00010',
                                               'helix-turn-helix transcriptional regulator_2', 'restriction_endonuclease', 
                                               'DNDJGHEP_00011', 'DNDJGHEP_00012', 'xerD', 'xerD/parA_A', 'xerD/parA_B',
                                               'parA', 'parB', 'pld', 'DNDJGHEP_00013', 'DNDJGHEP_00014', 'DNDJGHEP_00015',
                                               'DNDJGHEP_00016', 'DNDJGHEP_00017', 'radC', 'korC', 'DNDJGHEP_00018 (½)',
                                               'IS1-family_2', 'DNDJGHEP_00019', 'DNDJGHEP_00018 (2/2)', 'DNDJGHEP_00018 (2/2)/ccgA1_A',
                                               'DNDJGHEP_00018 (2/2)/ccgA1_B', 'ccgA1', 'DNDJGHEP_00020',
                                               'DNDJGHEP_00020/DNDJGHEP_00021_A', 'DNDJGHEP_00020/DNDJGHEP_00021_B',
                                               'DNDJGHEP_00021', 'ymoA', 'DNDJGHEP_00022', 'DNDJGHEP_00023', 'DNDJGHEP_00024',
                                               'DNDJGHEP_00025', 'klcA', 'DNDJGHEP_00026', 'DNDJGHEP_00027', 'ssb', 'DNDJGHEP_00028',
                                               'DNDJGHEP_00029', 'DNDJGHEP_00030', 'DNDJGHEP_00030/mobC_A',
                                               'DNDJGHEP_00030/mobC_B', 'mobC', 'mobC/nikA_A', 'mobC/nikA_B',
                                               'nikA', 'nikB', 'traH', 'traI', 'traJ', 'traK', 'traC', 'traL',
                                               'DNDJGHEP_00031', 'H-NS-like family', 'traM', 'traN', 'traO', 'traP',
                                               'traQ', 'traR', 'DNDJGHEP_00032', 'traU', 'traW',
                                               'traX-like', 'traY-like', 'excA', 'repC', 'repC/repA_A',
                                               'repC/repA_B')), y = factor(strain, levels = rev(c("K153", "PF_KPN15",
                                                                                              "PF_KPN18", "K163",
                                                                                              "PF_KPN11", "PF_KPN01",
                                                                                              "PF_KPN07", "PF_KPN10",
                                                                                              "PF_EC08", "PF_EC22",
                                                                                              "C288", "PF_EC23",
                                                                                              "J53", "PF_EC05"))), fill = median_log2FC)) +
  geom_tile() +
  facet_wrap(~treatment, nrow = 2) +
  xlab("pOXA-48 genes") +
  ylab("Strain") +
  scale_fill_gradient(low = "blue", high = "yellow", aesthetics = "fill") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))#, vjust = 0.5))

mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  ggplot(aes(x = annot_gene, y = strain, fill = median_log2FC)) +
  geom_tile() +
  xlab("pOXA-48") +
  ylab("Strain") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))

heatmap_mlog2FC <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  ggplot(aes(x = annot_gene, y = strain, fill = median_log2FC)) +
  facet_wrap(~treatment, nrow = 2) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "yellow") +
  theme_bw() +
  ylab("Strain") +
  xlab("pOXA-48 gene") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank())

heatmap_mlog2FC <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  ggplot(aes(x = annot_gene, y = strain, fill = median_log2FC)) +
  facet_wrap(~treatment, nrow = 2) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "yellow") +
  theme_bw() +
  ylab("Strain") +
  xlab("pOXA-48 gene") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(strip.background = element_blank())

ggplotly(heatmap_mlog2FC)

# Get only timepoint 2

mlog2FC <- mlog2FC %>%
  filter(timepoint == "t2")

mlog2FC %>%
  pivot_wider(names_from = annot_gene, values_from = median_log2FC)

mlog2FC %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_noab <- mlog2FC %>%
  filter(treatment == "LB+DAPG") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_ab <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

mlog2FC %>%
  ggplot(aes(x = factor(annot_gene, levels = c('repA', 'repA/trbC_A', 'repA/trbC_B',
                                               'trbC', 'trbB', 'trbA', 'trbN', 'DNDJGHEP_00001',
                                               'tir (1/2)', 'IS4-family_1', 'IS1-family_1', 'DNDJGHEP_00002',
                                               'DNDJGHEP_00002/blaOXA-48_A', 'DNDJGHEP_00002/blaOXA-48_B',
                                               'blaOXA-48', 'blaOXA-48/lysR_A', 'blaOXA-48/lysR_B', 'lysR',
                                               'lysR/IS4-family_A', 'lysR/IS4-family_B', 'IS4-family_2', 'tir (2/2)',
                                               'pemI', 'pemK', 'ltrA', 'ltrA/umuD_A', 'ltrA/umuD_B', 'umuD', 'umuC', 
                                               'DNDJGHEP_00003', 'DNDJGHEP_00004', 'DNDJGHEP_00005',
                                               'helix-turn-helix transcriptional regulator-1', 'relB',
                                               'relB/DNDJGHEP_00006_A', 'relB/DNDJGHEP_00006_B', 'DNDJGHEP_00006',
                                               'DNDJGHEP_00007', 'DNDJGHEP_00008', 'DNDJGHEP_00008/DNDJGHEP_00009_A',
                                               'DNDJGHEP_00008/DNDJGHEP_00009_B', 'DNDJGHEP_00009', 'DNDJGHEP_00010',
                                               'helix-turn-helix transcriptional regulator_2', 'restriction_endonuclease', 
                                               'DNDJGHEP_00011', 'DNDJGHEP_00012', 'xerD', 'xerD/parA_A', 'xerD/parA_B',
                                               'parA', 'parB', 'pld', 'DNDJGHEP_00013', 'DNDJGHEP_00014', 'DNDJGHEP_00015',
                                               'DNDJGHEP_00016', 'DNDJGHEP_00017', 'radC', 'korC', 'DNDJGHEP_00018 (½)',
                                               'IS1-family_2', 'DNDJGHEP_00019', 'DNDJGHEP_00018 (2/2)', 'DNDJGHEP_00018 (2/2)/ccgA1_A',
                                               'DNDJGHEP_00018 (2/2)/ccgA1_B', 'ccgA1', 'DNDJGHEP_00020',
                                               'DNDJGHEP_00020/DNDJGHEP_00021_A', 'DNDJGHEP_00020/DNDJGHEP_00021_B',
                                               'DNDJGHEP_00021', 'ymoA', 'DNDJGHEP_00022', 'DNDJGHEP_00023', 'DNDJGHEP_00024',
                                               'DNDJGHEP_00025', 'klcA', 'DNDJGHEP_00026', 'DNDJGHEP_00027', 'ssb', 'DNDJGHEP_00028',
                                               'DNDJGHEP_00029', 'DNDJGHEP_00030', 'DNDJGHEP_00030/mobC_A',
                                               'DNDJGHEP_00030/mobC_B', 'mobC', 'mobC/nikA_A', 'mobC/nikA_B',
                                               'nikA', 'nikB', 'traH', 'traI', 'traJ', 'traK', 'traC', 'traL',
                                               'DNDJGHEP_00031', 'H-NS-like family', 'traM', 'traN', 'traO', 'traP',
                                               'traQ', 'traR', 'DNDJGHEP_00032', 'traU', 'traW',
                                               'traX-like', 'traY-like', 'excA', 'repC', 'repC/repA_A',
                                               'repC/repA_B')), y = factor(strain, levels = rev(c("K153", "PF_KPN15",
                                                                                                  "PF_KPN18", "K163",
                                                                                                  "PF_KPN11", "PF_KPN01",
                                                                                                  "PF_KPN07", "PF_KPN10",
                                                                                                  "PF_EC08", "PF_EC22",
                                                                                                  "C288", "PF_EC23",
                                                                                                  "J53", "PF_EC05"))), fill = median_log2FC)) +
  geom_tile() +
  facet_wrap(~treatment, nrow = 2) +
  xlab("pOXA-48 genes") +
  ylab("Strain") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))#, vjust = 0.5))


# Get only timepoint 1

mlog2FC <- mlog2FC %>%
  filter(timepoint == "t1")

mlog2FC %>%
  pivot_wider(names_from = annot_gene, values_from = median_log2FC)

mlog2FC %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_noab <- mlog2FC %>%
  filter(treatment == "LB+DAPG") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

df_heatmap_ab <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  select(strain, annot_gene, median_log2FC) %>%
  reshape(idvar = "strain", timevar = "annot_gene", direction = "wide")

mlog2FC %>%
  filter(strain != "J53")%>%
  ggplot(aes(x = factor(annot_gene, levels = c('repA', 'repA/trbC_A', 'repA/trbC_B',
                                               'trbC', 'trbB', 'trbA', 'trbN', 'DNDJGHEP_00001',
                                               'tir (1/2)', 'IS4-family_1', 'IS1-family_1', 'DNDJGHEP_00002',
                                               'DNDJGHEP_00002/blaOXA-48_A', 'DNDJGHEP_00002/blaOXA-48_B',
                                               'blaOXA-48', 'blaOXA-48/lysR_A', 'blaOXA-48/lysR_B', 'lysR',
                                               'lysR/IS4-family_A', 'lysR/IS4-family_B', 'IS4-family_2', 'tir (2/2)',
                                               'pemI', 'pemK', 'ltrA', 'ltrA/umuD_A', 'ltrA/umuD_B', 'umuD', 'umuC', 
                                               'DNDJGHEP_00003', 'DNDJGHEP_00004', 'DNDJGHEP_00005',
                                               'helix-turn-helix transcriptional regulator-1', 'relB',
                                               'relB/DNDJGHEP_00006_A', 'relB/DNDJGHEP_00006_B', 'DNDJGHEP_00006',
                                               'DNDJGHEP_00007', 'DNDJGHEP_00008', 'DNDJGHEP_00008/DNDJGHEP_00009_A',
                                               'DNDJGHEP_00008/DNDJGHEP_00009_B', 'DNDJGHEP_00009', 'DNDJGHEP_00010',
                                               'helix-turn-helix transcriptional regulator_2', 'restriction_endonuclease', 
                                               'DNDJGHEP_00011', 'DNDJGHEP_00012', 'xerD', 'xerD/parA_A', 'xerD/parA_B',
                                               'parA', 'parB', 'pld', 'DNDJGHEP_00013', 'DNDJGHEP_00014', 'DNDJGHEP_00015',
                                               'DNDJGHEP_00016', 'DNDJGHEP_00017', 'radC', 'korC', 'DNDJGHEP_00018 (½)',
                                               'IS1-family_2', 'DNDJGHEP_00019', 'DNDJGHEP_00018 (2/2)', 'DNDJGHEP_00018 (2/2)/ccgA1_A',
                                               'DNDJGHEP_00018 (2/2)/ccgA1_B', 'ccgA1', 'DNDJGHEP_00020',
                                               'DNDJGHEP_00020/DNDJGHEP_00021_A', 'DNDJGHEP_00020/DNDJGHEP_00021_B',
                                               'DNDJGHEP_00021', 'ymoA', 'DNDJGHEP_00022', 'DNDJGHEP_00023', 'DNDJGHEP_00024',
                                               'DNDJGHEP_00025', 'klcA', 'DNDJGHEP_00026', 'DNDJGHEP_00027', 'ssb', 'DNDJGHEP_00028',
                                               'DNDJGHEP_00029', 'DNDJGHEP_00030', 'DNDJGHEP_00030/mobC_A',
                                               'DNDJGHEP_00030/mobC_B', 'mobC', 'mobC/nikA_A', 'mobC/nikA_B',
                                               'nikA', 'nikB', 'traH', 'traI', 'traJ', 'traK', 'traC', 'traL',
                                               'DNDJGHEP_00031', 'H-NS-like family', 'traM', 'traN', 'traO', 'traP',
                                               'traQ', 'traR', 'DNDJGHEP_00032', 'traU', 'traW',
                                               'traX-like', 'traY-like', 'excA', 'repC', 'repC/repA_A',
                                               'repC/repA_B')), y = factor(strain, levels = rev(c("K153", "PF_KPN15",
                                                                                                  "PF_KPN18", "K163",
                                                                                                  "PF_KPN11", "PF_KPN01",
                                                                                                  "PF_KPN07", "PF_KPN10",
                                                                                                  "PF_EC08", "PF_EC22",
                                                                                                  "C288", "PF_EC23",
                                                                                                  "J53", "PF_EC05"))), fill = median_log2FC)) +
  geom_tile() +
  facet_wrap(~treatment, nrow = 2) +
  xlab("pOXA-48 genes") +
  ylab("Strain") +
  scale_fill_gradient(low = "blue", high = "yellow") +
  theme(strip.background = element_blank()) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1))#, vjust = 0.5))
