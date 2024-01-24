#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 16:28:47 2023

@author: Jorge Sastre Domínguez

Script based on the code of F. Rousset available at https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb
for analysing the CRISPRi screenings based on gRNA counts per gene obtained with the CRISPR_merge_replicates.py script.

This script takes the guides count table obtained by CRISPR_merge_replicates.py and performs the pooling of the data
and calculous of Log2FC and gene scores

"""

import pandas as pd
from collections import Counter
import itertools
from scipy.stats.stats import pearsonr,spearmanr
import os
from tqdm import tqdm
import re
#from adjustText import *
from Bio import SeqIO
import scipy.stats as st
from Bio.Seq import Seq
from scipy.cluster.hierarchy import dendrogram, linkage  
from Bio.SeqUtils import nt_search
import sklearn.metrics as metrics
#import statsmodels.api as sm
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)

# Import library  and get control guide sequences

library = pd.read_excel('/home/jorge/Documents/CRISPRi/library_controls_final.xlsx')
control_guides=list(library[library.name.str.match('control')]['guide'])

# Import data and metadata
data = pd.read_excel("/home/jorge/Documents/CRISPRi/guides_count/merged_counts_per_sample.xlsx", sheet_name= "counts per sample")
pooled_data = pd.read_excel("/home/jorge/Documents/CRISPRi/guides_count/merged_counts_per_sample.xlsx", sheet_name= "counts per sample")

# We import again data and metadata code to have the metadata table from the raw count guide table per replicate
# And our pooled data is the table merged by sample obtained by the previous script CRISPR_merge_replicates.py

# Create metadata table from the sample names of the merged counts

metadata = pd.DataFrame(columns=["sample", "strain", "timepoint", "replicate", "treatment"])

merged_counts = pd.read_excel("/home/jorge/Documents/CRISPRi/guides_count/merged_counts.xlsx", sheet_name= "merged counts")

samplenames = list(merged_counts.columns[1:])
samples = []
strainnames = []
timepoints = []
replicates = []
treatments = []
screens = []

for sample in samplenames[:]:
    strainname = "_".join(sample.split("_")[:-2])
    strainnames.append(strainname)
    
    timepoint = sample.split("_")[-2]
    timepoints.append(timepoint)
    
    replicate = sample.split("_")[-1]
    replicates.append(replicate)
    
    if sample[0] == "e":
        treatment = "LB+DAPG+ERTA"
        treatments.append(treatment)
    else:
        treatment = "LB+DAPG"
        treatments.append(treatment)
    
    screen = treatment+"."+timepoint
    screens.append(screen)
    
    samples.append(strainname+"_"+timepoint)

metadata = pd.DataFrame(list(zip(samples, samplenames, strainnames, timepoints, replicates, treatments, screens)),
                        columns = ["sample_ID", "name", "strain", "timepoint", "replicate", "treatment", "screen"])

#print(metadata)

# Set guides as row indexes as done in the original code for both data and pooled data dfs

#print(pooled_data)
data.set_index("Guide", drop = True, inplace = True)
#print(pooled_data)
data.index.names = [None]

#print(pooled_data)
pooled_data.set_index("Guide", drop = True, inplace = True)
#print(pooled_data)
pooled_data.index.names = [None]

#Remove guides where all T0 samples have less than 20 reads
pooled_data=pooled_data[pooled_data.filter(regex='_t0').max(axis=1)>20]

#For each strain, set guides that have <20 reads in T0 sample as NaN.
samples=[c for c in pooled_data.columns if ('_t0' not in c)]
strains=list(set([s for s in samples]))

for s in strains:
    if s[0] != "e":
        pooled_data.loc[data[s[:-1]+'0']<20,[c for c in pooled_data.columns if (s in c)&('_t0' not in c)]]=np.nan
    else:
        pooled_data.loc[data[s[1:-1]+'0']<20,[c for c in pooled_data.columns if (s in c)&('_t0' not in c)]]=np.nan

# Normalize read counts by sample size
norm_factors=np.mean(pooled_data.mean())/pooled_data.mean()
pooled_data=pooled_data*norm_factors
#print(norm_factors)

# In our case, we don't need to check for RM systems, and perfect matches have already been analysed during library construction

#### NOW WE CALCULATE LOG2FC VALUES
#print(metadata)

log2FC=pd.DataFrame(index=pooled_data.index)
strains = metadata.sample_ID.unique()
for strain in strains:
    if strain[0] == "e": # Loop for ERTA condition
        for screen in [c for c in pooled_data.columns if (strain == c)&('_t0' not in c)]: #Select strain samples (except t0) and calculate log2FC for each
            #print(screen+ " ERTA")    
            FC=pd.DataFrame({screen:np.log2((pooled_data[screen]+1)/(pooled_data[strain[1:-1]+'0']+1))},index=pooled_data.index)
            log2FC=pd.concat([log2FC,FC],axis=1,sort=False)
    else: # Loop for No ERTA condition
        for screen in [c for c in pooled_data.columns if (strain == c)&('_t0' not in c)]: #Select strain samples (except t0) and calculate log2FC for each
            #print(screen + " No ERTA")    
            FC=pd.DataFrame({screen:np.log2((pooled_data[screen]+1)/(pooled_data[strain[:-1]+'0']+1))},index=pooled_data.index)
            log2FC=pd.concat([log2FC,FC],axis=1,sort=False)  
#Center log2FC values by removing the median of control guides
controls_median=log2FC.loc[control_guides].median()
log2FC=log2FC-controls_median

# Now we need to merge the guides of the Log2FC table with the corresponding genes of the library
# Set guide as column to merge
log2FC['guide']=log2FC.index
log2FC=log2FC.loc[:,sorted(log2FC.columns)]
first_column = log2FC.pop("guide")
log2FC.insert(0, "guide", first_column)
log2FC=pd.merge(log2FC,library,on='guide')

# Add gene column by removing suffix of library names (gene_1, gene_2, etc.)
log2FC["gene"] = log2FC["name"].str[:-2]

# Now, lets calculate median log2FC by gene

medians=log2FC[samples+['gene']].groupby('gene').median()

# Save medians table into xlsx for R representation
"""
outname = "/home/jorge/Documents/CRISPRi/median_log2FC_by_gene.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
medians.to_excel(writer, sheet_name='median log2FC by gene', index=True)
metadata.to_excel(writer, sheet_name="metadata", index = False)
writer.save()
"""
# Save scores per guide to check which works better table into xlsx for R representation
"""
outname = "/home/jorge/Documents/CRISPRi/median_log2FC_by_guide.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
log2FC.to_excel(writer, sheet_name='median log2FC by gene', index=True)
#metadata.to_excel(writer, sheet_name="metadata", index = False)
writer.save()
"""
"""
# Create scatterplots for all samples
# For each strain we'll include the 6 timepoints, 3 per treatment
for s in metadata.strain.unique():
    #print(s)
    strain_plot = []
    if s[0]!="e":
        for tested_s in [c for c in medians.columns if s in c]:
            strain_plot.append(tested_s)
        #print(strain_plot)
        df_strain_plot = medians[medians.columns.intersection(strain_plot)]
        df_strain_plot_t = df_strain_plot.T
        
        plt.figure(figsize = (7,4))
        sns.scatterplot(data = df_strain_plot_t, s = 100, legend = None)
        plt.axhline(y=0)
        plt.xlabel("Timepoint & Condition")
        plt.ylabel("Median gene score")
        plt.xticks(rotation=45, ha = "right")
        plt.title(s)
        plt.show()
#plt.savefig()"""