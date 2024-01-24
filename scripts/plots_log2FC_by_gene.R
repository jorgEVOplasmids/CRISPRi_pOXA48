
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

#median_log2FC <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene.xlsx", 1)
#metadata <- read.xlsx("/home/jorge/Documents/CRISPRi/results/median_log2FC_by_gene.xlsx", 2)
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

#highlight_blaoxa <- mlog2FC_w_metadata %>%
  #filter(strain == "C288" | strain == "eC288") %>%
  #filter(gene == "blaOXA-48")

mypal <- RColorBrewer::brewer.pal(7, "Set1")
mypal <- mypal[-c(5,6,7)]
mypal <- append(mypal, c("#808080","#ffdf22"))

scores_by_strain <- mlog2FC %>%
  #filter(strain == "C288" | strain == "eC288") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.1)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  facet_wrap(~strain, nrow = 2)+
  xlab("Treatment + Timepoint")+
  ylab("Gene score (Log2 FC)")+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(scores_by_strain)

C288_scores <- mlog2FC %>%
  filter(strain == "C288") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("C288")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(C288_scores)

J53_scores <- mlog2FC %>%
  filter(strain == "J53") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("J53")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(J53_scores)

K153_scores <- mlog2FC %>%
  filter(strain == "K153") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("K153")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(K153_scores)

K163_scores <- mlog2FC %>%
  filter(strain == "K163") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("K163")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(K163_scores)

PF_EC05_scores <- mlog2FC %>%
  filter(strain == "PF_EC05") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_EC05")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_EC05_scores)

PF_EC08_scores <- mlog2FC %>%
  filter(strain == "PF_EC05") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_EC08")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_EC08_scores)

PF_EC22_scores <- mlog2FC %>%
  filter(strain == "PF_EC22") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_EC22")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_EC22_scores)

PF_EC23_scores <- mlog2FC %>%
  filter(strain == "PF_EC23") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_EC23")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_EC23_scores)

PF_KPN01_scores <- mlog2FC %>%
  filter(strain == "PF_KPN01") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN01")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN01_scores)

PF_KPN07_scores <- mlog2FC %>%
  filter(strain == "PF_KPN07") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN07")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN07_scores)

PF_KPN10_scores <- mlog2FC %>%
  filter(strain == "PF_KPN10") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN10")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN10_scores)

PF_KPN11_scores <- mlog2FC %>%
  filter(strain == "PF_KPN11") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN11")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN11_scores)

PF_KPN15_scores <- mlog2FC %>%
  filter(strain == "PF_KPN15") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN15")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN15_scores)

PF_KPN18_scores <- mlog2FC %>%
  filter(strain == "PF_KPN18") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = annot_gene, col = annotation)) +
  geom_jitter(width = 0.25)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Timepoint + Treatment")+
  ylab("Gene score (Log2 FC)")+
  ggtitle("PF_KPN18")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 22)+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(PF_KPN18_scores)

#### NOW, REPRESENT FOR EACH STRAIN AND TIMEPOINT, ERTA VS NO ERTA

# Change sample names in table to remove the e at the beginning which indicates the treatment already shown in the treatment column

for(i in 1:nrow(`mlog2FC`)) {       # for-loop over rows
  if (str_sub(mlog2FC$sample_ID[i], start = 1, end = length(mlog2FC$sample_ID[i])) == "e") {
    mlog2FC$sample_ID[i] <- str_sub(mlog2FC$sample_ID[i], start = 2)
  }
}

# So we get the scores with antibiotics in another column for plotting
erta_mlog2FC <- mlog2FC %>% filter(treatment == "LB+DAPG+ERTA") %>% arrange(sample_ID)
erta_mlog2FC$genebysample <- interaction(erta_mlog2FC$annot_gene, erta_mlog2FC$sample_ID)
no_erta_mlog2FC <- mlog2FC %>% filter(treatment == "LB+DAPG") %>% arrange(sample_ID)
no_erta_mlog2FC$genebysample <- interaction(no_erta_mlog2FC$annot_gene, no_erta_mlog2FC$sample_ID)
erta_mlog2FC <- erta_mlog2FC %>% select(median_log2FC, genebysample)
no_erta_mlog2FC <- merge(no_erta_mlog2FC, erta_mlog2FC, by = "genebysample")

#no_erta_mlog2FC$erta_median_log2FC <- erta_mlog2FC$median_log2FC

colnames(no_erta_mlog2FC)[5] = "median_log2FC"
colnames(no_erta_mlog2FC)[11] = "erta_median_log2FC"
no_erta_mlog2FC = subset(no_erta_mlog2FC, select = -c(1) )


# All strains

strains_erta_noerta <- no_erta_mlog2FC %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("ERTA vs NO ERTA")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 12)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(strains_erta_noerta)

# And plot for each strain

erta_vs_no_erta_C288 <- no_erta_mlog2FC %>%
  filter(strain == "C288") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("C288")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_C288)

erta_vs_no_erta_J53 <- no_erta_mlog2FC %>%
  filter(strain == "J53") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("J53")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_J53)

erta_vs_no_erta_K153 <- no_erta_mlog2FC %>%
  filter(strain == "K153") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("K153")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_K153)

erta_vs_no_erta_K163 <- no_erta_mlog2FC %>%
  filter(strain == "K163") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("K163")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_K163)

erta_vs_no_erta_PF_EC05 <- no_erta_mlog2FC %>%
  filter(strain == "PF_EC05") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_EC05")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_EC05)

erta_vs_no_erta_PF_EC08 <- no_erta_mlog2FC %>%
  filter(strain == "PF_EC08") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_EC08")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_EC08)

erta_vs_no_erta_PF_EC22 <- no_erta_mlog2FC %>%
  filter(strain == "PF_EC22") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_EC22")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_EC22)

erta_vs_no_erta_PF_EC23 <- no_erta_mlog2FC %>%
  filter(strain == "PF_EC23") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_EC23")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_EC23)

erta_vs_no_erta_PF_KPN01 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN01") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN01")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN01)

erta_vs_no_erta_PF_KPN07 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN07") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN07")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN07)

erta_vs_no_erta_PF_KPN10 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN10") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN10")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN10)

erta_vs_no_erta_PF_KPN11 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN11") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN11")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN11)

erta_vs_no_erta_PF_KPN15 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN15") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN15")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN15)

erta_vs_no_erta_PF_KPN18 <- no_erta_mlog2FC %>%
  filter(strain == "PF_KPN18") %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = annot_gene, col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  ggtitle("PF_KPN18")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~strain+timepoint, ncol = 6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(erta_vs_no_erta_PF_KPN18)

# Save table with information of erta and no_erta

#write_xlsx(no_erta_mlog2FC, "~/Documents/CRISPRi/table_erta_vs_no_erta.xlsx")

# Plot by species to check differences

species_results <- mlog2FC %>%
  #filter(strain == "C288" | strain == "eC288") %>%
  ggplot(aes(x = interaction(timepoint, treatment), y = median_log2FC, group = interaction(annot_gene, strain), col = annotation)) +
  geom_jitter(width = 0.1)+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  facet_wrap(~species, nrow = 2)+
  xlab("Treatment + Timepoint")+
  ylab("Gene score (Log2 FC)")+
  ylim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(species_results)

species_results_erta_vs_no_erta <- no_erta_mlog2FC %>%
  ggplot(aes(x = median_log2FC, y = erta_median_log2FC, group = interaction(annot_gene, strain), col = annotation)) +
  geom_point()+
  #geom_jitter(data = highlight_blaoxa, aes(x = sample_ID, y = median_log2FC, group = sample_ID), col = "red", size = 2) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
  geom_vline(xintercept = 0, linetype = 2, col = "darkgrey") +
  xlab("Gene score (LB+DAPG)")+
  ylab("Gene score (LB+DAPG+ERTA)")+
  #geom_text(data = filter(gene == "blaOXA-48"), aes(label = gene))+
  theme_bw(base_size = 16)+
  ylim(-10.5,10.5)+
  xlim(-10.5,10.5)+
  scale_color_manual(values = mypal)+
  facet_wrap(~species+timepoint, nrow = 2)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

ggplotly(species_results_erta_vs_no_erta)

