
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
library(lmPerm)
library(mev)

#median_log2FC <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene.xlsx", 1)
#metadata <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene.xlsx", 2)
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

### Perform permutation test to check significant log2FC values in the distribution of our gene scores

# Plot histogram of the distribution of the gene scores at timepoint t3

mlog2FC %>%
  filter(timepoint == "t3") %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# Plot histogram of the distribution of the gene scores at timepoint t3 and in presence of ERTA

mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA", timepoint == "t3") %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# Split df by strain and plot distribution of scores
#####
# C288

mlog2FC_C288 <- mlog2FC %>%
  filter(strain == "C288", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_C288_ERTA <- mlog2FC %>%
  filter(strain == "C288", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_C288 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_C288_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# J53

mlog2FC_J53 <- mlog2FC %>%
  filter(strain == "J53", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_J53_ERTA <- mlog2FC %>%
  filter(strain == "J53", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_J53 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_J53_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# K153

mlog2FC_K153 <- mlog2FC %>%
  filter(strain == "K153", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_K153_ERTA <- mlog2FC %>%
  filter(strain == "K153", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_K153 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_K153_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# K163

mlog2FC_K163 <- mlog2FC %>%
  filter(strain == "K163", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_K163_ERTA <- mlog2FC %>%
  filter(strain == "K163", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_K163 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_K163_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_EC05

mlog2FC_PF_EC05 <- mlog2FC %>%
  filter(strain == "PF_EC05", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_EC05_ERTA <- mlog2FC %>%
  filter(strain == "PF_EC05", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_EC05 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_EC05_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_EC08

mlog2FC_PF_EC08 <- mlog2FC %>%
  filter(strain == "PF_EC08", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_EC08_ERTA <- mlog2FC %>%
  filter(strain == "PF_EC08", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_EC08 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_EC08_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_EC22

mlog2FC_PF_EC22 <- mlog2FC %>%
  filter(strain == "PF_EC22", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_EC22_ERTA <- mlog2FC %>%
  filter(strain == "PF_EC22", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_EC22 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_EC22_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_EC23

mlog2FC_PF_EC23 <- mlog2FC %>%
  filter(strain == "PF_EC23", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_EC23_ERTA <- mlog2FC %>%
  filter(strain == "PF_EC23", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_EC23 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_EC23_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN01

mlog2FC_PF_KPN01 <- mlog2FC %>%
  filter(strain == "PF_KPN01", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN01_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN01", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN01 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN01_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN07

mlog2FC_PF_KPN07 <- mlog2FC %>%
  filter(strain == "PF_KPN07", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN07_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN07", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN07 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN07_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN10

mlog2FC_PF_KPN10 <- mlog2FC %>%
  filter(strain == "PF_KPN10", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN10_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN10", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN10 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN10_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN11

mlog2FC_PF_KPN11 <- mlog2FC %>%
  filter(strain == "PF_KPN11", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN11_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN11", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN11 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN11_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN15

mlog2FC_PF_KPN15 <- mlog2FC %>%
  filter(strain == "PF_KPN15", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN15_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN15", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN15 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN15_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

# PF_KPN18

mlog2FC_PF_KPN18 <- mlog2FC %>%
  filter(strain == "PF_KPN18", timepoint == "t3", treatment == "LB+DAPG")

mlog2FC_PF_KPN18_ERTA <- mlog2FC %>%
  filter(strain == "PF_KPN18", timepoint == "t3", treatment == "LB+DAPG+ERTA")

mlog2FC_PF_KPN18 %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

mlog2FC_PF_KPN18_ERTA %>%
  ggplot(aes(x = median_log2FC)) +
  geom_histogram(position = "identity")

#####

# Check if each gene score is significantly different from the rest of the scores in each condition

# With all the strains together

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results$treatment <- "LB+DAPG"
stat_results$p_adj <- p.adjust(stat_results$p_valors, method = "fdr")


#### Now with antibiotic

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  p_value_gene <- sum(perm_stat > obs_stat)/num_permutations
  pvals <- c(pvals, p_value_gene)
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_ERTA <- as.data.frame(cbind(unique(table$annot_gene), p_valors))
stat_results_ERTA$treatment <- "LB+DAPG+ERTA"
stat_results_ERTA$p_adj <- p.adjust(stat_results_ERTA$p_valors, method = "fdr")

############################

# Save tables of stats results as a xlsx file

write_xlsx(stat_results, "permutation_test_merged_no_ERTA.xlsx")
write_xlsx(stat_results_ERTA, "permutation_test_merged_ERTA.xlsx")

############################

# For each of the species

# E. coli without ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", species == "E. coli")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results$treatment <- "LB+DAPG"
stat_results$p_adj <- p.adjust(stat_results$p_valors, method = "fdr")

#### E. coli with ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", species == "E. coli")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  p_value_gene <- sum(perm_stat > obs_stat)/num_permutations
  pvals <- c(pvals, p_value_gene)
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_ERTA <- as.data.frame(cbind(unique(table$annot_gene), p_valors))
stat_results_ERTA$treatment <- "LB+DAPG+ERTA"
stat_results_ERTA$p_adj <- p.adjust(stat_results_ERTA$p_valors, method = "fdr")

# Save tables of stats results as a xlsx file

write_xlsx(stat_results, "permutation_test_merged_no_ERTA_ecoli.xlsx")
write_xlsx(stat_results_ERTA, "permutation_test_merged_ERTA_ecoli.xlsx")

# K. pneumoniae without ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", species == "K. pneumoniae")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results$treatment <- "LB+DAPG"
stat_results$p_adj <- p.adjust(stat_results$p_valors, method = "fdr")

#### K. pneumoniae with ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", species == "K. pneumoniae")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(mean(shuffled_group[1:num_observations]) - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  p_value_gene <- sum(perm_stat > obs_stat)/num_permutations
  pvals <- c(pvals, p_value_gene)
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_ERTA <- as.data.frame(cbind(unique(table$annot_gene), p_valors))
stat_results_ERTA$treatment <- "LB+DAPG+ERTA"
stat_results_ERTA$p_adj <- p.adjust(stat_results_ERTA$p_valors, method = "fdr")

# Save tables of stats results as a xlsx file

write_xlsx(stat_results, "permutation_test_merged_no_ERTA_kleb.xlsx")
write_xlsx(stat_results_ERTA, "permutation_test_merged_ERTA_kleb.xlsx")

############################

# For each strain individually

# C288
#####
table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "C288")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_C288 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_C288$p_valors <- as.numeric(stat_results_C288$p_valors)

#stat_results_C288$padj <- p.adjust(stat_results_C288$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "C288")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_C288_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_C288_ERTA$p_valors <- as.numeric(stat_results_C288_ERTA$p_valors)

#stat_results_C288_ERTA$padj <- p.adjust(stat_results_C288_ERTA$p_valors, method = "fdr")
#####
# J53

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "J53")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_J53 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_J53$p_valors <- as.numeric(stat_results_J53$p_valors)

#stat_results_J53$padj <- p.adjust(stat_results_J53$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "J53")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_J53_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_J53_ERTA$p_valors <- as.numeric(stat_results_J53_ERTA$p_valors)

#stat_results_J53_ERTA$padj <- p.adjust(stat_results_J53_ERTA$p_valors, method = "fdr")
#####

#####
# K153

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "K153")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_K153 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_K153$p_valors <- as.numeric(stat_results_K153$p_valors)

#stat_results_K153$padj <- p.adjust(stat_results_K153$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "K153")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_K153_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_K153_ERTA$p_valors <- as.numeric(stat_results_K153_ERTA$p_valors)

#stat_results_K153_ERTA$padj <- p.adjust(stat_results_K153_ERTA$p_valors, method = "fdr")
#####

#####
# K163

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "K163")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_K163 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_K163$p_valors <- as.numeric(stat_results_K163$p_valors)

#stat_results_K163$padj <- p.adjust(stat_results_K163$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "K163")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_K163_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_K163_ERTA$p_valors <- as.numeric(stat_results_K163_ERTA$p_valors)

#stat_results_K163_ERTA$padj <- p.adjust(stat_results_K163_ERTA$p_valors, method = "fdr")

#####

#####
# PF_KPN01

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN01")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN01 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN01$p_valors <- as.numeric(stat_results_PF_KPN01$p_valors)

#stat_results_PF_KPN01$padj <- p.adjust(stat_results_PF_KPN01$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN01")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN01_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN01_ERTA$p_valors <- as.numeric(stat_results_PF_KPN01_ERTA$p_valors)

#stat_results_PF_KPN01_ERTA$padj <- p.adjust(stat_results_PF_KPN01_ERTA$p_valors, method = "fdr")
#####

#####
# PF_KPN07

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN07")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN07 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN07$p_valors <- as.numeric(stat_results_PF_KPN07$p_valors)

#stat_results_PF_KPN07$padj <- p.adjust(stat_results_PF_KPN07$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN07")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN07_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN07_ERTA$p_valors <- as.numeric(stat_results_PF_KPN07_ERTA$p_valors)

#stat_results_PF_KPN01_ERTA$padj <- p.adjust(stat_results_PF_KPN01_ERTA$p_valors, method = "fdr")

#####

#####
# PF_KPN10

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN10")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN10 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN10$p_valors <- as.numeric(stat_results_PF_KPN10$p_valors)

#stat_results_PF_KPN10$padj <- p.adjust(stat_results_PF_KPN10$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN10")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN10_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN10_ERTA$p_valors <- as.numeric(stat_results_PF_KPN10_ERTA$p_valors)

#stat_results_PF_KPN01_ERTA$padj <- p.adjust(stat_results_PF_KPN01_ERTA$p_valors, method = "fdr")

#####

# PF_KPN11

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN11")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN11 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN11$p_valors <- as.numeric(stat_results_PF_KPN11$p_valors)

#stat_results_PF_KPN11$padj <- p.adjust(stat_results_PF_KPN11$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN11")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN11_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN11_ERTA$p_valors <- as.numeric(stat_results_PF_KPN11_ERTA$p_valors)

#stat_results_PF_KPN11_ERTA$padj <- p.adjust(stat_results_PF_KPN11_ERTA$p_valors, method = "fdr")

#####

# PF_KPN15

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN15")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN15 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN15$p_valors <- as.numeric(stat_results_PF_KPN15$p_valors)

#stat_results_PF_KPN15$padj <- p.adjust(stat_results_PF_KPN15$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN15")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN15_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN15_ERTA$p_valors <- as.numeric(stat_results_PF_KPN15_ERTA$p_valors)

#stat_results_PF_KPN15_ERTA$padj <- p.adjust(stat_results_PF_KPN15_ERTA$p_valors, method = "fdr")

#####

#####

# PF_KPN18

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_KPN18")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN18 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN18$p_valors <- as.numeric(stat_results_PF_KPN18$p_valors)

#stat_results_PF_KPN18$padj <- p.adjust(stat_results_PF_KPN18$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_KPN18")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_KPN18_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_KPN18_ERTA$p_valors <- as.numeric(stat_results_PF_KPN18_ERTA$p_valors)

#stat_results_PF_KPN18_ERTA$padj <- p.adjust(stat_results_PF_KPN18_ERTA$p_valors, method = "fdr")

#####

#####

# PF_EC05

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_EC05")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC05 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC05$p_valors <- as.numeric(stat_results_PF_EC05$p_valors)

#stat_results_PF_EC05$padj <- p.adjust(stat_results_PF_EC05$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_EC05")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC05_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC05_ERTA$p_valors <- as.numeric(stat_results_PF_EC05_ERTA$p_valors)

#stat_results_PF_EC05_ERTA$padj <- p.adjust(stat_results_PF_EC05_ERTA$p_valors, method = "fdr")

#####

# PF_EC08

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_EC08")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC08 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC08$p_valors <- as.numeric(stat_results_PF_EC08$p_valors)

#stat_results_PF_EC08$padj <- p.adjust(stat_results_PF_EC08$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_EC08")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC08_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC08_ERTA$p_valors <- as.numeric(stat_results_PF_EC08_ERTA$p_valors)

#stat_results_PF_EC08_ERTA$padj <- p.adjust(stat_results_PF_EC08_ERTA$p_valors, method = "fdr")


# PF_EC22

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_EC22")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC22 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC22$p_valors <- as.numeric(stat_results_PF_EC22$p_valors)

#stat_results_PF_EC22$padj <- p.adjust(stat_results_PF_EC22$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_EC22")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC22_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC22_ERTA$p_valors <- as.numeric(stat_results_PF_EC22_ERTA$p_valors)

#stat_results_PF_EC22_ERTA$padj <- p.adjust(stat_results_PF_EC22_ERTA$p_valors, method = "fdr")

# PF_EC23

table <- mlog2FC %>% filter(treatment == "LB+DAPG", timepoint == "t3", strain == "PF_EC23")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC23 <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC23$p_valors <- as.numeric(stat_results_PF_EC23$p_valors)

#stat_results_PF_EC23$padj <- p.adjust(stat_results_PF_EC23$p_valors, method = "fdr")

# With ERTA

table <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", timepoint == "t3", strain == "PF_EC23")
pvals <- c()
p_valors <- c()

for (g in unique(table$annot_gene)) {
  #print(g)
  
  # Get diference between gene score and the mean of the rest of the scores as the observed statistic
  
  gene_score <- subset(table, annot_gene == g)$median_log2FC
  rest_scores <- subset(table, annot_gene != g)$median_log2FC
  
  obs_stat <- abs(gene_score - mean(rest_scores))
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 100000
  num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    shuffled_group <- sample(c(gene_score, rest_scores))
    perm_stat[i, ] <- abs(shuffled_group[1:num_observations] - mean(shuffled_group[(num_observations + 1):length(shuffled_group)]))
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the score of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values)
  
}

stat_results_PF_EC23_ERTA <- as.data.frame(cbind(unique(table$annot_gene),p_valors))
stat_results_PF_EC23_ERTA$p_valors <- as.numeric(stat_results_PF_EC23_ERTA$p_valors)

#stat_results_PF_EC23_ERTA$padj <- p.adjust(stat_results_PF_EC23_ERTA$p_valors, method = "fdr")

##### Represent Volcano plots for stats results

# Build table for representation. For global stats, calculate mean score per gene (as we will show individual points, one per gene)
# and merge with the table with adjusted p values for each of the conditions. Finally, represent.

# No ERTA

mlog_2FC_no_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG")
mean_scores_no_ERTA <- aggregate(mlog_2FC_no_ERTA$median_log2FC, by = list(mlog_2FC_no_ERTA$annot_gene), FUN = mean)
colnames(mean_scores_no_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_no_ERTA <- merge(stat_results, mean_scores_no_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test

table_volcano_no_ERTA[table_volcano_no_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_no_ERTA, lab = table_volcano_no_ERTA$annot_gene,
                x = "log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG')

# ERTA

mlog_2FC_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA")
mean_scores_ERTA <- aggregate(mlog_2FC_ERTA$median_log2FC, by = list(mlog_2FC_ERTA$annot_gene), FUN = mean)
colnames(mean_scores_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_ERTA) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_ERTA <- merge(stat_results_ERTA, mean_scores_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test

table_volcano_ERTA[table_volcano_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_ERTA, lab = table_volcano_ERTA$annot_gene,
                x = "log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA')

#####

# For each strain

## C288

mlog_2FC_C288 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "C288", timepoint == "t3")
#mean_scores_C288 <- aggregate(mlog_2FC_C288$median_log2FC, by = list(mlog_2FC_C288$annot_gene), FUN = mean)
#colnames(mean_scores_C288) <- c("annot_gene", "log2FC")
colnames(stat_results_C288) <- c("annot_gene", "p_valors")
table_volcano_C288 <- merge(stat_results, mlog_2FC_C288, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_C288$p_valors <- as.numeric(table_volcano_C288$p_valors)
table_volcano_C288[table_volcano_C288$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_C288, lab = table_volcano_C288$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (C288)')

# ERTA

mlog_2FC_C288_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "C288", timepoint == "t3")
#mean_scores_C288_ERTA <- aggregate(mlog_2FC_C288_ERTA$median_log2FC, by = list(mlog_2FC_C288_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_C288_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_C288_ERTA) <- c("annot_gene", "p_valors")
table_volcano_C288_ERTA <- merge(stat_results_C288_ERTA, mlog_2FC_C288_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_C288_ERTA$p_valors <- as.numeric(table_volcano_C288_ERTA$p_valors)
table_volcano_C288_ERTA[table_volcano_C288_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_C288_ERTA, lab = table_volcano_C288_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (C288)')

## J53

mlog_2FC_J53 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "J53", timepoint == "t3")
#mean_scores_J53 <- aggregate(mlog_2FC_J53$median_log2FC, by = list(mlog_2FC_J53$annot_gene), FUN = mean)
#colnames(mean_scores_J53) <- c("annot_gene", "log2FC")
colnames(stat_results_J53) <- c("annot_gene", "p_valors")
table_volcano_J53 <- merge(stat_results, mlog_2FC_J53, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_J53$p_valors <- as.numeric(table_volcano_J53$p_valors)
table_volcano_J53[table_volcano_J53$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_J53, lab = table_volcano_J53$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (J53)')

# ERTA

mlog_2FC_J53_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "J53", timepoint == "t3")
#mean_scores_J53_ERTA <- aggregate(mlog_2FC_J53_ERTA$median_log2FC, by = list(mlog_2FC_J53_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_J53_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_J53_ERTA) <- c("annot_gene", "p_valors")
table_volcano_J53_ERTA <- merge(stat_results_J53_ERTA, mlog_2FC_J53_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_J53_ERTA$p_valors <- as.numeric(table_volcano_J53_ERTA$p_valors)
table_volcano_J53_ERTA[table_volcano_J53_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_J53_ERTA, lab = table_volcano_J53_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (J53)')

## K153

mlog_2FC_K153 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "K153", timepoint == "t3")
#mean_scores_K153 <- aggregate(mlog_2FC_K153$median_log2FC, by = list(mlog_2FC_K153$annot_gene), FUN = mean)
#colnames(mean_scores_K153) <- c("annot_gene", "log2FC")
colnames(stat_results_K153) <- c("annot_gene", "p_valors")
table_volcano_K153 <- merge(stat_results, mlog_2FC_K153, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_K153$p_valors <- as.numeric(table_volcano_K153$p_valors)
table_volcano_K153[table_volcano_K153$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_K153, lab = table_volcano_K153$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (K153)')

# ERTA

mlog_2FC_K153_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "K153", timepoint == "t3")
#mean_scores_K153_ERTA <- aggregate(mlog_2FC_K153_ERTA$median_log2FC, by = list(mlog_2FC_K153_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_K153_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_K153_ERTA) <- c("annot_gene", "p_valors")
table_volcano_K153_ERTA <- merge(stat_results_K153_ERTA, mlog_2FC_K153_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_K153_ERTA$p_valors <- as.numeric(table_volcano_K153_ERTA$p_valors)
table_volcano_K153_ERTA[table_volcano_K153_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_K153_ERTA, lab = table_volcano_K153_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (K153)')

## K163

mlog_2FC_K163 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "K163", timepoint == "t3")
#mean_scores_K163 <- aggregate(mlog_2FC_K163$median_log2FC, by = list(mlog_2FC_K163$annot_gene), FUN = mean)
#colnames(mean_scores_K163) <- c("annot_gene", "log2FC")
colnames(stat_results_K163) <- c("annot_gene", "p_valors")
table_volcano_K163 <- merge(stat_results, mlog_2FC_K163, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_K163$p_valors <- as.numeric(table_volcano_K163$p_valors)
table_volcano_K163[table_volcano_K163$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_K163, lab = table_volcano_K163$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (K163)')

# ERTA

mlog_2FC_K163_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "K163", timepoint == "t3")
#mean_scores_K163_ERTA <- aggregate(mlog_2FC_K163_ERTA$median_log2FC, by = list(mlog_2FC_K163_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_K163_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_K163_ERTA) <- c("annot_gene", "p_valors")
table_volcano_K163_ERTA <- merge(stat_results_K163_ERTA, mlog_2FC_K163_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_K163_ERTA$p_valors <- as.numeric(table_volcano_K163_ERTA$p_valors)
table_volcano_K163_ERTA[table_volcano_K163_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_K163_ERTA, lab = table_volcano_K163_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (K163)')

## PF_KPN01

mlog_2FC_PF_KPN01 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN01", timepoint == "t3")
#mean_scores_PF_KPN01 <- aggregate(mlog_2FC_PF_KPN01$median_log2FC, by = list(mlog_2FC_PF_KPN01$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN01) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN01) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN01 <- merge(stat_results, mlog_2FC_PF_KPN01, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN01$p_valors <- as.numeric(table_volcano_PF_KPN01$p_valors)
table_volcano_PF_KPN01[table_volcano_PF_KPN01$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN01, lab = table_volcano_PF_KPN01$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN01)')

# ERTA

mlog_2FC_PF_KPN01_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN01", timepoint == "t3")
#mean_scores_PF_KPN01_ERTA <- aggregate(mlog_2FC_PF_KPN01_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN01_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN01_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN01_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN01_ERTA <- merge(stat_results_PF_KPN01_ERTA, mlog_2FC_PF_KPN01_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN01_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN01_ERTA$p_valors)
table_volcano_PF_KPN01_ERTA[table_volcano_PF_KPN01_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN01_ERTA, lab = table_volcano_PF_KPN01_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN01)')

## PF_KPN07

mlog_2FC_PF_KPN07 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN07", timepoint == "t3")
#mean_scores_PF_KPN07 <- aggregate(mlog_2FC_PF_KPN07$median_log2FC, by = list(mlog_2FC_PF_KPN07$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN07) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN07) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN07 <- merge(stat_results, mlog_2FC_PF_KPN07, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN07$p_valors <- as.numeric(table_volcano_PF_KPN07$p_valors)
table_volcano_PF_KPN07[table_volcano_PF_KPN07$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN07, lab = table_volcano_PF_KPN07$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN07)')

# ERTA

mlog_2FC_PF_KPN07_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN07", timepoint == "t3")
#mean_scores_PF_KPN07_ERTA <- aggregate(mlog_2FC_PF_KPN07_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN07_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN07_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN07_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN07_ERTA <- merge(stat_results_PF_KPN07_ERTA, mlog_2FC_PF_KPN07_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN07_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN07_ERTA$p_valors)
table_volcano_PF_KPN07_ERTA[table_volcano_PF_KPN07_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN07_ERTA, lab = table_volcano_PF_KPN07_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN07)')

## PF_KPN10

mlog_2FC_PF_KPN10 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN10", timepoint == "t3")
#mean_scores_PF_KPN10 <- aggregate(mlog_2FC_PF_KPN10$median_log2FC, by = list(mlog_2FC_PF_KPN10$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN10) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN10) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN10 <- merge(stat_results, mlog_2FC_PF_KPN10, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN10$p_valors <- as.numeric(table_volcano_PF_KPN10$p_valors)
table_volcano_PF_KPN10[table_volcano_PF_KPN10$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN10, lab = table_volcano_PF_KPN10$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN10)')

# ERTA

mlog_2FC_PF_KPN10_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN10", timepoint == "t3")
#mean_scores_PF_KPN10_ERTA <- aggregate(mlog_2FC_PF_KPN10_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN10_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN10_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN10_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN10_ERTA <- merge(stat_results_PF_KPN10_ERTA, mlog_2FC_PF_KPN10_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN10_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN10_ERTA$p_valors)
table_volcano_PF_KPN10_ERTA[table_volcano_PF_KPN10_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN10_ERTA, lab = table_volcano_PF_KPN10_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN10)')

## PF_KPN11

mlog_2FC_PF_KPN11 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN11", timepoint == "t3")
#mean_scores_PF_KPN11 <- aggregate(mlog_2FC_PF_KPN11$median_log2FC, by = list(mlog_2FC_PF_KPN11$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN11) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN11) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN11 <- merge(stat_results, mlog_2FC_PF_KPN11, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN11$p_valors <- as.numeric(table_volcano_PF_KPN11$p_valors)
table_volcano_PF_KPN11[table_volcano_PF_KPN11$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN11, lab = table_volcano_PF_KPN11$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN11)')

# ERTA

mlog_2FC_PF_KPN11_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN11", timepoint == "t3")
#mean_scores_PF_KPN11_ERTA <- aggregate(mlog_2FC_PF_KPN11_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN11_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN11_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN11_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN11_ERTA <- merge(stat_results_PF_KPN11_ERTA, mlog_2FC_PF_KPN11_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN11_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN11_ERTA$p_valors)
table_volcano_PF_KPN11_ERTA[table_volcano_PF_KPN11_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN11_ERTA, lab = table_volcano_PF_KPN11_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN11)')

## PF_KPN15

mlog_2FC_PF_KPN15 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN15", timepoint == "t3")
#mean_scores_PF_KPN15 <- aggregate(mlog_2FC_PF_KPN15$median_log2FC, by = list(mlog_2FC_PF_KPN15$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN15) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN15) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN15 <- merge(stat_results, mlog_2FC_PF_KPN15, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN15$p_valors <- as.numeric(table_volcano_PF_KPN15$p_valors)
table_volcano_PF_KPN15[table_volcano_PF_KPN15$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN15, lab = table_volcano_PF_KPN15$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN15)')

# ERTA

mlog_2FC_PF_KPN15_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN15", timepoint == "t3")
#mean_scores_PF_KPN15_ERTA <- aggregate(mlog_2FC_PF_KPN15_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN15_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN15_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN15_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN15_ERTA <- merge(stat_results_PF_KPN15_ERTA, mlog_2FC_PF_KPN15_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN15_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN15_ERTA$p_valors)
table_volcano_PF_KPN15_ERTA[table_volcano_PF_KPN15_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN15_ERTA, lab = table_volcano_PF_KPN15_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN15)')

## PF_KPN18

mlog_2FC_PF_KPN18 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_KPN18", timepoint == "t3")
#mean_scores_PF_KPN18 <- aggregate(mlog_2FC_PF_KPN18$median_log2FC, by = list(mlog_2FC_PF_KPN18$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN18) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN18) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN18 <- merge(stat_results, mlog_2FC_PF_KPN18, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN18$p_valors <- as.numeric(table_volcano_PF_KPN18$p_valors)
table_volcano_PF_KPN18[table_volcano_PF_KPN18$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN18, lab = table_volcano_PF_KPN18$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_KPN18)')

# ERTA

mlog_2FC_PF_KPN18_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_KPN18", timepoint == "t3")
#mean_scores_PF_KPN18_ERTA <- aggregate(mlog_2FC_PF_KPN18_ERTA$median_log2FC, by = list(mlog_2FC_PF_KPN18_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_KPN18_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_KPN18_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_KPN18_ERTA <- merge(stat_results_PF_KPN18_ERTA, mlog_2FC_PF_KPN18_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_KPN18_ERTA$p_valors <- as.numeric(table_volcano_PF_KPN18_ERTA$p_valors)
table_volcano_PF_KPN18_ERTA[table_volcano_PF_KPN18_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_KPN18_ERTA, lab = table_volcano_PF_KPN18_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_KPN18)')

## PF_EC05

mlog_2FC_PF_EC05 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_EC05", timepoint == "t3")
#mean_scores_PF_EC05 <- aggregate(mlog_2FC_PF_EC05$median_log2FC, by = list(mlog_2FC_PF_EC05$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC05) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC05) <- c("annot_gene", "p_valors")
table_volcano_PF_EC05 <- merge(stat_results, mlog_2FC_PF_EC05, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC05$p_valors <- as.numeric(table_volcano_PF_EC05$p_valors)
table_volcano_PF_EC05[table_volcano_PF_EC05$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_EC05, lab = table_volcano_PF_EC05$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_EC05)')

# ERTA

mlog_2FC_PF_EC05_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_EC05", timepoint == "t3")
#mean_scores_PF_EC05_ERTA <- aggregate(mlog_2FC_PF_EC05_ERTA$median_log2FC, by = list(mlog_2FC_PF_EC05_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC05_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC05_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_EC05_ERTA <- merge(stat_results_PF_EC05_ERTA, mlog_2FC_PF_EC05_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC05_ERTA$p_valors <- as.numeric(table_volcano_PF_EC05_ERTA$p_valors)
table_volcano_PF_EC05_ERTA[table_volcano_PF_EC05_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_EC05_ERTA, lab = table_volcano_PF_EC05_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_EC05)')

## PF_EC08

mlog_2FC_PF_EC08 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_EC08", timepoint == "t3")
#mean_scores_PF_EC08 <- aggregate(mlog_2FC_PF_EC08$median_log2FC, by = list(mlog_2FC_PF_EC08$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC08) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC08) <- c("annot_gene", "p_valors")
table_volcano_PF_EC08 <- merge(stat_results, mlog_2FC_PF_EC08, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC08$p_valors <- as.numeric(table_volcano_PF_EC08$p_valors)
table_volcano_PF_EC08[table_volcano_PF_EC08$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_EC08, lab = table_volcano_PF_EC08$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_EC08)')

# ERTA

mlog_2FC_PF_EC08_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_EC08", timepoint == "t3")
#mean_scores_PF_EC08_ERTA <- aggregate(mlog_2FC_PF_EC08_ERTA$median_log2FC, by = list(mlog_2FC_PF_EC08_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC08_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC08_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_EC08_ERTA <- merge(stat_results_PF_EC08_ERTA, mlog_2FC_PF_EC08_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC08_ERTA$p_valors <- as.numeric(table_volcano_PF_EC08_ERTA$p_valors)
table_volcano_PF_EC08_ERTA[table_volcano_PF_EC08_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_EC08_ERTA, lab = table_volcano_PF_EC08_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_EC08)')

## PF_EC22

mlog_2FC_PF_EC22 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_EC22", timepoint == "t3")
#mean_scores_PF_EC22 <- aggregate(mlog_2FC_PF_EC22$median_log2FC, by = list(mlog_2FC_PF_EC22$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC22) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC22) <- c("annot_gene", "p_valors")
table_volcano_PF_EC22 <- merge(stat_results, mlog_2FC_PF_EC22, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC22$p_valors <- as.numeric(table_volcano_PF_EC22$p_valors)
table_volcano_PF_EC22[table_volcano_PF_EC22$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_EC22, lab = table_volcano_PF_EC22$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_EC22)')

# ERTA

mlog_2FC_PF_EC22_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_EC22", timepoint == "t3")
#mean_scores_PF_EC22_ERTA <- aggregate(mlog_2FC_PF_EC22_ERTA$median_log2FC, by = list(mlog_2FC_PF_EC22_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC22_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC22_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_EC22_ERTA <- merge(stat_results_PF_EC22_ERTA, mlog_2FC_PF_EC22_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC22_ERTA$p_valors <- as.numeric(table_volcano_PF_EC22_ERTA$p_valors)
table_volcano_PF_EC22_ERTA[table_volcano_PF_EC22_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_EC22_ERTA, lab = table_volcano_PF_EC22_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_EC22)')

## PF_EC23

mlog_2FC_PF_EC23 <- mlog2FC %>% filter(treatment == "LB+DAPG", strain == "PF_EC23", timepoint == "t3")
#mean_scores_PF_EC23 <- aggregate(mlog_2FC_PF_EC23$median_log2FC, by = list(mlog_2FC_PF_EC23$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC23) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC23) <- c("annot_gene", "p_valors")
table_volcano_PF_EC23 <- merge(stat_results_PF_EC23, mlog_2FC_PF_EC23, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC23$p_valors <- as.numeric(table_volcano_PF_EC23$p_valors)
table_volcano_PF_EC23[table_volcano_PF_EC23$p_valors == 0, ]$p_valors <- 1/100000

EnhancedVolcano(table_volcano_PF_EC23, lab = table_volcano_PF_EC23$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (PF_EC23)')

# ERTA

mlog_2FC_PF_EC23_ERTA <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", strain == "PF_EC23", timepoint == "t3")
#mean_scores_PF_EC23_ERTA <- aggregate(mlog_2FC_PF_EC23_ERTA$median_log2FC, by = list(mlog_2FC_PF_EC23_ERTA$annot_gene), FUN = mean)
#colnames(mean_scores_PF_EC23_ERTA) <- c("annot_gene", "log2FC")
colnames(stat_results_PF_EC23_ERTA) <- c("annot_gene", "p_valors")
table_volcano_PF_EC23_ERTA <- merge(stat_results_PF_EC23_ERTA, mlog_2FC_PF_EC23_ERTA, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_PF_EC23_ERTA$p_valors <- as.numeric(table_volcano_PF_EC23_ERTA$p_valors)
table_volcano_PF_EC23_ERTA[table_volcano_PF_EC23_ERTA$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_PF_EC23_ERTA, lab = table_volcano_PF_EC23_ERTA$annot_gene,
                x = "median_log2FC", y = "p_valors", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (PF_EC23)')

# Volcano plots for data of species

# E. coli

# No ERTA

stat_results_ecoli <- read.xlsx("permutation_test_merged_no_ERTA_ecoli.xlsx", sheetIndex = 1)

mlog_2FC_ecoli <- mlog2FC %>% filter(treatment == "LB+DAPG", species == "E. coli", timepoint == "t3")
mean_scores_no_ERTA_ecoli <- aggregate(mlog_2FC_ecoli$median_log2FC, by = list(mlog_2FC_ecoli$annot_gene), FUN = mean)
colnames(mean_scores_no_ERTA_ecoli) <- c("annot_gene", "median_log2FC")
colnames(stat_results_ecoli) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_no_ERTA_ecoli <- merge(stat_results_ecoli, mean_scores_no_ERTA_ecoli, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_no_ERTA_ecoli$p_valors <- as.numeric(table_volcano_no_ERTA_ecoli$p_valors)
table_volcano_no_ERTA_ecoli[table_volcano_no_ERTA_ecoli$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_no_ERTA_ecoli, lab = table_volcano_no_ERTA_ecoli$annot_gene,
                x = "median_log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (E. coli)')



# ERTA

stat_results_ecoli <- read.xlsx("permutation_test_merged_ERTA_ecoli.xlsx", sheetIndex = 1)

mlog_2FC_ecoli <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", species == "E. coli", timepoint == "t3")
mean_scores_ERTA_ecoli <- aggregate(mlog_2FC_ecoli$median_log2FC, by = list(mlog_2FC_ecoli$annot_gene), FUN = mean)
colnames(mean_scores_ERTA_ecoli) <- c("annot_gene", "median_log2FC")
colnames(stat_results_ecoli) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_ERTA_ecoli <- merge(stat_results_ecoli, mean_scores_ERTA_ecoli, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_ERTA_ecoli$p_valors <- as.numeric(table_volcano_ERTA_ecoli$p_valors)
table_volcano_ERTA_ecoli[table_volcano_ERTA_ecoli$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_ERTA_ecoli, lab = table_volcano_ERTA_ecoli$annot_gene,
                x = "median_log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (E. coli)')



# K. pneumoniae

# No ERTA

stat_results_kleb <- read.xlsx("permutation_test_merged_no_ERTA_kleb.xlsx", sheetIndex = 1)

mlog_2FC_kleb <- mlog2FC %>% filter(treatment == "LB+DAPG", species == "K. pneumoniae", timepoint == "t3")
mean_scores_no_ERTA_kleb <- aggregate(mlog_2FC_kleb$median_log2FC, by = list(mlog_2FC_kleb$annot_gene), FUN = mean)
colnames(mean_scores_no_ERTA_kleb) <- c("annot_gene", "median_log2FC")
colnames(stat_results_kleb) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_no_ERTA_kleb <- merge(stat_results_kleb, mean_scores_no_ERTA_kleb, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_no_ERTA_kleb$p_valors <- as.numeric(table_volcano_no_ERTA_kleb$p_valors)
table_volcano_no_ERTA_kleb[table_volcano_no_ERTA_kleb$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_no_ERTA_kleb, lab = table_volcano_no_ERTA_kleb$annot_gene,
                x = "median_log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG (K. pneumoniae)')



# ERTA

stat_results_kleb <- read.xlsx("permutation_test_merged_ERTA_kleb.xlsx", sheetIndex = 1)

mlog_2FC_kleb <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA", species == "K. pneumoniae", timepoint == "t3")
mean_scores_ERTA_kleb <- aggregate(mlog_2FC_kleb$median_log2FC, by = list(mlog_2FC_kleb$annot_gene), FUN = mean)
colnames(mean_scores_ERTA_kleb) <- c("annot_gene", "median_log2FC")
colnames(stat_results_kleb) <- c("annot_gene", "p_valors", "treatment", "p_adj")
table_volcano_ERTA_kleb <- merge(stat_results_kleb, mean_scores_ERTA_kleb, by = "annot_gene")

# replace p values equal to 0 with the minimum p value we can get by the permutation test
table_volcano_ERTA_kleb$p_valors <- as.numeric(table_volcano_ERTA_kleb$p_valors)
table_volcano_ERTA_kleb[table_volcano_ERTA_kleb$p_adj == 0, ]$p_adj <- 1/100000

EnhancedVolcano(table_volcano_ERTA_kleb, lab = table_volcano_ERTA_kleb$annot_gene,
                x = "median_log2FC", y = "p_adj", pCutoff = 0.05, FCcutoff = 0.5, labSize = 3,
                colAlpha = 4/5,
                legendPosition = 'top',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black',
                title = 'LB+DAPG+ERTA (K. pneumoniae)')

