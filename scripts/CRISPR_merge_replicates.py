#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:20:04 2023

@author: Jorge Sastre Domínguez

Script based on the code of F. Rousset available at https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb
for analysing the CRISPRi screenings based on gRNA counts per gene obtained with the merge_counts.py script.

This first script basically takes the merged table of the guide counts, and checks the reproducibility between
replicates. After plotting the Pearson's R coefficient, it merges all replicates into their same sample
to perform subsequent analyses. 

"""

import pandas as pd
from scipy.stats.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import statistics

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

# Rename counts file

data = merged_counts

##### CALCULATE EXPERIMENTAL REPRODUCIBILITY

strains=metadata.strain.unique()
corr_coefs=pd.DataFrame(columns=['strain','screen','coef'])

for condition in ["LB+DAPG", "LB+DAPG+ERTA"]:
    strains=list(metadata.loc[metadata.treatment==condition,'sample_ID'].values) # Keep strains assayed in this condition
    #print(strains)
    strains=list(set([s for s in strains if strains.count(s)==3])) # Keep samples with 3 replicates
    for s in strains:
        strain_samples = [s+"_r1",s+"_r2", s+"_r3"] # Get the 3 replicates names
        #print(strain_samples)
        corr_coefs.loc[len(corr_coefs)]=[s,condition,pearsonr(data[strain_samples[0]],data[strain_samples[1]])[0]] #Compute Pearson's R
        
#print(corr_coefs)

for i,screen in enumerate(["LB+DAPG", "LB+DAPG+ERTA"]):
    fig=plt.figure(figsize=(4,1.5))
    ax = fig.add_subplot(111)
    plt.hist(corr_coefs.loc[corr_coefs.screen==screen,'coef'],color=sns.color_palette()[i],edgecolor='white',bins=np.arange(0.9,1.1,0.01))
    plt.xlim([0.8,1])
    plt.xlabel('Pearson\'s R of biological replicates')
    plt.ylabel('Sample triplicates')
    plt.yticks([0,5,10])
    plt.text(x=0.05,y=0.9,s=screen,transform=ax.transAxes,color=sns.color_palette()[i])
    sns.despine()
    plt.show()
    image_format = 'svg' # e.g .png, .svg, etc.
    image_name = 'pearson_corr_'+str(screen)+'.svg'
    fig.savefig(image_name, format=image_format, dpi=1200)
    
print('Median correlation coefficient : '+str(round(corr_coefs.coef.median(),3)))

# Now that we have checked the reproducibility of the biological triplicates, we can merge the guides_count xlsx file columns
# into sample columns rather than replicate columns. Then, the counts for each sample will be the sum of the gRNAs counted
# in each replicate. Then, coverage will be calculated based on these values for each sample

#print(metadata.loc[metadata.treatment=="LB+DAPG+ERTA"].sample.unique())

# List of samples IDs without replicates, only strain+timepoint

#print(len(metadata.sample_ID.unique()))

sampleIDlist = metadata.sample_ID.unique() # Kind of the same phylosophy as before with each strain
data_by_sample = {} # Dictionary which will keep each sample column as key (sample name): value (sum of 3 replicates guide counts)
data_by_sample["Guide"] = data["Guide"] # Get guides columns and include as first column of new dataframe

for sID in sampleIDlist: # Initiate loop for each sample (strain+timepoint)
    #print(sID)
    sID_replicates = []
    for column in data.columns: # Iterate through the data columns (merged_counts.csv)
        if column == sID+"_r1" or column == sID+"_r2" or column == sID+"_r3":
            #print(sID, column)
            #sID_replicates.append(column)
            sID_replicates.append(column) # In sID_replicates list we include the replicates of each sample
    data_by_sample[sID] = data[sID_replicates].sum(axis = 1)
    #print(sID, data[sID_replicates].sum(axis = 1))
    #print(sID, sID_replicates)

merged_counts_by_sample = pd.DataFrame.from_dict(data_by_sample, orient = "columns")

# Now, get the sheet with the same info but including gene names and coverage by sample

merged_counts_gene = pd.read_excel("/home/jorge/Documents/CRISPRi/guides_count/merged_counts.xlsx", sheet_name= "merged counts with gene names")
merged_counts_gene_by_sample = pd.DataFrame.from_dict(data_by_sample, orient = "columns")
merged_counts_gene_by_sample["name"] = merged_counts_gene["name"]

# Move gene name to first position

first_column = merged_counts_gene_by_sample.pop("name")
merged_counts_gene_by_sample.insert(0, "name", first_column)

cov_row = ["", "Sample coverage per gene"]
# Declare length of guides library to calculate coverage as the median number of counts per guide multiplied by the number of guides per gene
lenlibrary = len(merged_counts_gene_by_sample.index)
for column in merged_counts_gene_by_sample:
    if column != "Guide" and column != "name":
        mediancountguides = statistics.median(merged_counts_gene_by_sample[column].values)
        covsample = (mediancountguides*5)
        #print(covsample)
        cov_row.append(covsample)

# Include coverage row in the second sheet
merged_counts_gene_by_sample.loc[len(merged_counts_gene_by_sample)] = cov_row
"""
counts_directory = "/home/jorge/Documents/CRISPRi/guides_count/"
# Save in an excel df with two identical sheets, the second one containing the gene names, as in the previous step but with the replicates merged
print("Saving xlsx tables")
outname = counts_directory+"merged_counts_per_sample.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
merged_counts_by_sample.to_excel(writer, sheet_name='counts per sample', index=False)
merged_counts_gene_by_sample.to_excel(writer, sheet_name='counts per sample gene names', index=False)
writer.save()

print("Output saved to "+counts_directory+outname)
"""