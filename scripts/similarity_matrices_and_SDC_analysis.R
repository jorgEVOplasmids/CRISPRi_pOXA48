
setwd("/home/jorge/Documents/CRISPRi/correlations_strains")

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
library(ggplot2)
library(reshape2)
library(corrplot)
library(ggcorrplot)

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

cost_info <- read.xlsx("/home/jorge/Documents/CRISPRi/correlations_costs/costes_curvas_cepas.xlsx", 1)

mlog2FC <- merge(mlog2FC, cost_info, by = "strain")
mlog2FC$coste_OXA <- as.numeric(mlog2FC$coste_OXA)
mlog2FC$median_log2FC <- as.numeric(mlog2FC$median_log2FC)

# Get only timepoint 3

mlog2FC <- mlog2FC %>%
  filter(timepoint == "t3")

# Compare whether the response is conserved among strains by performing a chi squared test
# Reshape the table to get a matrix

table_int_noERTA <- mlog2FC %>%
  filter(treatment == "LB+DAPG") %>%
  select(strain, annot_gene, median_log2FC)
table_int_noERTA$median_log2FC <- as.numeric(table_int_noERTA$median_log2FC)

# Obtain matrix for statistics

table_int_noERTA <- as.data.frame(table_int_noERTA)
table_int_noERTA <- table_int_noERTA %>% pivot_wider(names_from = strain, values_from = median_log2FC)
table_int_noERTA <- as.data.frame(table_int_noERTA)
rownames(table_int_noERTA) <- table_int_noERTA$annot_gene
table_int_noERTA <- table_int_noERTA %>% select(-c(annot_gene))

cor_matrix <- cor(table_int_noERTA)

corplot_noERTA <- corrplot(cor_matrix, type = "lower", method = "circle", tl.col = "black")

write.xlsx(cor_matrix, "correlation_matrix_noERTA.xlsx")

####### Same with ERTA

# Reshape the table to get a matrix

table_int_noERTA <- mlog2FC %>%
  filter(treatment == "LB+DAPG+ERTA") %>%
  select(strain, annot_gene, median_log2FC)
table_int_noERTA$median_log2FC <- as.numeric(table_int_noERTA$median_log2FC)
#colnames(table_int_noERTA) <- c("strain", "annot_gene", "median_log2FC_noERTA")

# Obtain matrix for statistics

table_int_noERTA <- as.data.frame(table_int_noERTA)
table_int_noERTA <- table_int_noERTA %>% pivot_wider(names_from = strain, values_from = median_log2FC)
table_int_noERTA <- as.data.frame(table_int_noERTA)
rownames(table_int_noERTA) <- table_int_noERTA$annot_gene
table_int_noERTA <- table_int_noERTA %>% select(-c(annot_gene))

cor_matrix <- cor(table_int_noERTA)

corrplot(cor_matrix, type = "lower", method = "circle", tl.col = "black")

write.xlsx(cor_matrix, "correlation_matrix_ERTA.xlsx")

### Code to calculate and plot SDC


# Function to calculate Sorensen-Dice similarity coefficient
sorensen_dice <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  return((2 * intersection) / (length(set1) + length(set2)))
}

# Calculate similarity coefficients
similarity_results <- df %>%
  group_by(Species, Condition) %>%
  summarise(similarity = list(combn(unique(Strain), 2, simplify = FALSE))) %>%
  unnest(similarity) %>%
  mutate(Strain1 = map_chr(similarity, 1),
         Strain2 = map_chr(similarity, 2)) %>%
  rowwise() %>%
  mutate(Similarity = sorensen_dice(
    df$Significant_genes[df$Strain == Strain1 & df$Condition == Condition & df$Species == Species],
    df$Significant_genes[df$Strain == Strain2 & df$Condition == Condition & df$Species == Species]
  )) %>%
  ungroup()

# Plot the similarity coefficients
ggplot(similarity_results, aes(x = Condition, y = Similarity)) +
  geom_boxplot(aes(fill = Species))+
  geom_jitter(aes(col = Species), width = 0.1) +
  ylim(0,1)+
  #facet_wrap(~ Condition + Species, scales = "free_x") +
  labs(x = "Strain Pair",
       y = "Sørensen-Dice Similarity Coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw(base_size = 20) +
  theme(panel.background = element_blank(), panel.grid = element_blank(),
        #strip.text.x = element_blank(),
        #legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,  hjust=1))

