# Initial set-up ----

setwd("/home/jorge/Documents/CRISPRi/delta_blaoxa")

library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggpmisc)
library(readr)
library(data.table)
library(purrr)
library(dplyr)
library(tidyr)
library(growthrates)
library (gggenes)
library(nlme)
library(ape)
library(ggrepel)
library(flux)
library(broom)
library(tibble)
library(lme4)
library(afex)
library(knitr)
library(ggpubr)
library(fmsb)
library(scales)
library(car)
library(rstatix)
library(mgcv)
library(performance)
library(see)
library(xlsx)

path_to_txt= "/home/jorge/Documents/CRISPRi/delta_blaoxa/txt_files/"

path_to_Growthrates= "/home/jorge/Documents/CRISPRi/delta_blaoxa/generated_files/"

path_to_output= "/home/jorge/Documents/CRISPRi/delta_blaoxa/generated_files/"

Plantilla <- read.table("/home/jorge/Documents/CRISPRi/delta_blaoxa/Plantilla.csv", na.strings="",header=TRUE, sep = ";")

# Preparing and organising the data ----
  
file.list <- list.files(path =path_to_txt, full.names = F)

df.list <- lapply(paste0(path_to_txt, file.list), 
                      function(x)read.delim(x, header=T, dec=","))

attr(df.list, "names") <- file.list

# Ahora le meto los tiempos (24 horas que son 1440 minutos)
# Cambio esta version un poco porque el PR2 no mide cada 10 min, sino cada 20...

t1 <- rep(seq(0, 1450-10, 10),1)
t2 <- rep(seq(0, 1450-10, 10),1)
t3 <- rep(seq(0, 1430-10, 10),1)
t4 <- rep(seq(0, 1420-10, 10),1)
t5 <- rep(seq(0, 1450-10, 10),1)
t6 <- rep(seq(0, 1450-10, 10),1)
t7 <- rep(seq(0, 1450-10, 10),1)
time <- c(t1,t2,t3,t4,t5,t6,t7)
df_data <- bind_rows(df.list, .id = "id") %>% 
  mutate(Time=time)

# Ahora le cambio el nombre de la optical density

df_test<-df_data %>% 
  select( -`TOpticalDensity600`) %>% 
  gather(-Time,-id,  key = Well, value = OD ) %>% 
  separate(id, into=c("Plate", "Date"), remove = F) 

# Ahora quito los "-" y lo presento a la Plantilla

curve_data <- df_test %>% 
  left_join(Plantilla %>% mutate(Plate=as.character(Plate), 
                                     Date=as.character(Date),
                                     Well=as.character(Well)))%>%
  mutate(Date=as.numeric(Date))%>%
  filter(Sample!="-")

# Data analysis ----

### Import new data

curve_data <- read.csv("curves_blaOXA.csv")

## AUC ----

# First I'll calculate the Area Under the Curves (AUC) with the flux package

data_analysed_AUC<- curve_data %>% 
  group_by(Plate , Sample ,  Clon , Replicate, Date,id,Plasmid, Strain, Treatment) %>% 
  group_modify(~ as.data.frame(flux::auc(.x$Time, .x$OD))) %>%
  mutate(AUC=`flux::auc(.x$Time, .x$OD)`) %>% 
  select(-`flux::auc(.x$Time, .x$OD)`)%>% 
  ungroup() 

# Calculate the mean AUC for the plasmid-free strain of each of the samples (takes into consideration sample and id, therefore if there is an effect of a day, it will take that into consideration (problem: wt and TC need to be in the same plate, same day))
# AUC_PF is the average of the AUC of the plasmid-free strain, organised by Sample
# AUC_Rel is the AUC of the TC normalised to the AUC_PF of each Sample, Repeat

info_data<-data_analysed_AUC %>% select(Sample,Replicate,Plasmid,Treatment,id,AUC,Strain,Date) %>% unique()

new_data<-info_data %>%
  filter(Plasmid=="none") %>%
  group_by(Sample,id,Plasmid) %>%
  mutate(AUC_PF = mean(AUC)) %>%
  ungroup()%>%
  select(AUC_PF,id,Strain,Date)

info_data2<-info_data %>% full_join(new_data) %>% unique()

data_analysed_AUC_RELATIVE<- info_data2 %>%
  mutate(AUCRel= AUC / AUC_PF) %>% unique()

data_analysed_AUC_RELATIVE %>% 
  ggplot(aes(Strain, AUC))+
  geom_boxplot(aes(x=Strain,color=Plasmid),notch=F,outlier.shape=T,show.legend = TRUE) +
  #geom_point(aes(alpha=0.3, color=Plasmid),shape=1)+
  #This takes the individual points out
  theme_linedraw()+
  labs(y=" AUC", x="Strain ",fill="Strain")+
  geom_hline(yintercept=1, linetype="dashed", 
             color = "grey", size=0.5)+
  ylim(700,1250)
  #coord_flip()

data_analysed_AUC_RELATIVE %>% 
  filter(Plasmid != "none") %>%
  ggplot(aes(Strain, AUCRel))+
  geom_boxplot(aes(x=Strain,color=Plasmid),notch=F,outlier.shape=T,show.legend = TRUE) +
  #geom_point(aes(alpha=0.3, color=Plasmid),shape=1)+
  #This takes the individual points out
  theme_linedraw()+
  labs(y=" AUC", x="Strain ",fill="Strain")+
  geom_hline(yintercept=1, linetype="dashed", 
             color = "grey", size=0.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())
#coord_flip()

# Merge all data with mean and standard error of each strain and condition

meantable <- data_analysed_AUC_RELATIVE %>%
  group_by(Sample, Plasmid) %>% 
  summarise(mean = mean(AUCRel),
            std = sd(AUCRel))

meantable <- as.data.frame(meantable)

# Extract order of samples by pOXA-48 fitness effect

poxamean <- meantable %>%
  filter(Plasmid == "pOXA-48")

orderedpoxamean <- poxamean[order(poxamean$mean),]

orderedpoxamean$Sample

# Represent barplot

meantable %>%
  filter(Plasmid != "none") %>%
  ggplot(aes(x = factor(Sample, levels = c(orderedpoxamean$Sample)), y = mean-1, fill = Plasmid, group = Plasmid))+
  #geom_boxplot(aes(x=Strain,color=Plasmid),notch=F,outlier.shape=T,show.legend = TRUE) +
  #geom_point(aes(alpha=0.3, color=Plasmid),shape=1)+
  #This takes the individual points out
  geom_bar(width=0.7, stat = "identity", position = position_dodge(width=0.8)) +
  theme_bw(base_size = 18)+
  ylim(-0.25, 0.25)+
  geom_errorbar(aes(x = factor(Sample, levels = c(orderedpoxamean$Sample)),ymax=mean-1+std, ymin = mean-1-std), position=position_dodge(width=0.8), colour="darkgrey", width = 0, size=0.7) +
  #facet_wrap(~Plasmid) +
  labs(y=" Relative AUC", x="Strain ",fill="Strain")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())

#write.xlsx(meantable, "costes_cepas_IC90.xlsx")

costes_IC90 <- read.xlsx("costes_cepas_IC90.xlsx", sheetIndex = 2)

costes_IC90 %>% 
  ggplot(aes(x = log(IC90), y = mean))+
  #geom_boxplot(aes(x=Strain,color=Plasmid),notch=F,outlier.shape=T,show.legend = TRUE) +
  geom_point()+
  stat_cor(method = "spearman", label.x.npc = 0.8) +
  geom_smooth(method = "loess", alpha = 0.1, col = "black") +  #This takes the individual points out
  #geom_bar(width=0.7, stat = "identity", position = position_dodge(width=0.8)) +
  theme_bw(base_size = 18)+
  #geom_errorbar(aes(x = factor(Sample, levels = c(orderedpoxamean$Sample)),ymax=mean-1+std, ymin = mean-1-std), position=position_dodge(width=0.8), colour="grey", width = 0, size=0.7) +
  #facet_wrap(~Plasmid) +
  labs(y=" Cost", x="Log(IC90) ")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.background = element_blank(), panel.grid = element_blank(), 
        #aspect.ratio = 1, 
        #legend.position = "none",
        strip.background = element_blank())


# Calculate statistical significance of fitness cost for each strain genotype vs their PF relative

table_pvalues <- data.frame(Strain=character(0), pvalue=numeric(0))

### pOXA-48 vs plasmid-free

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN11",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN11" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN11" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJN46",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJN46" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJN46" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN10",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN10" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN10" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K219",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K219" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K219" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC728",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC728" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC728" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "C288",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "C288" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "C288" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "EC22",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "EC22" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "EC22" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "CF12",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "CF12" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "CF12" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC310",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC310" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC310" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "N23",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "N23" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "N23" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN18",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN18" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN18" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC527",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC527" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC527" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K163",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K163" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K163" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC164",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC164" & data_analysed_AUC_RELATIVE$Plasmid == "pOXA-48", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC164" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)


### deltablaoxa-48_pOXA-48 vs plasmid-free

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN11d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN11" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN11" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJN46d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJN46" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJN46" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN10d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN10" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN10" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC[1:7],
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K219d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K219" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K219" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC728d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC728" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC728" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "C288d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "C288" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "C288" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "EC22d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "EC22" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "EC22" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "CF12d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "CF12" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "CF12" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC310d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC310" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC310" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "N23d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "N23" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "N23" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "KPN18d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN18" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "KPN18" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC527d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC527" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC527" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "K163d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K163" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "K163" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues <- table_pvalues %>% add_row (Strain= "AJC164d",
                                            pvalue= t.test(data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC164" & data_analysed_AUC_RELATIVE$Plasmid == "delta_blaOXA", ]$AUC,
                                                           data_analysed_AUC_RELATIVE[data_analysed_AUC_RELATIVE$Strain == "AJC164" & data_analysed_AUC_RELATIVE$Plasmid == "none", ]$AUC,
                                                           paired=TRUE)$p.value)

table_pvalues$p.adj.fdr <- p.adjust(table_pvalues$pvalue, method="bonferroni", n=28)

