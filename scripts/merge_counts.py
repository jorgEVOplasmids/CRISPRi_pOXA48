#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 16:42:53 2023

@author: Jorge Sastre Domínguez

Script for merging every CRISPRi count table from individual samples into the same table, being
each row a guide, each column a sample, and each cell the number of occurrences of the guide in the sample.

"""

import pandas as pd
import os
from functools import reduce
import statistics

counts_directory = "/home/jorge/Documents/CRISPRi/guides_count/"
counts_filenames = os.listdir(counts_directory)

tablelist = []

for i, element in enumerate(counts_filenames):
    
    filetable = pd.read_csv(counts_directory+counts_filenames[i], sep = " ", header=0)
    
    # For each guide in each sample, get the number of counts
    
    samplecounts = filetable.columns[0]
    
    # Remove guides with less counts than 5 at t0 and 2 at other times to avoid df size error when merging
    if "t0" in counts_filenames[i]:
        filetable_t0 = filetable[filetable[samplecounts] > 5]
        tablelist.append(filetable_t0)
    else:
        filetable_t = filetable[filetable[samplecounts] > 2]
        tablelist.append(filetable_t)
    
    print("Appending "+counts_filenames[i])
    
# Merge everything in the same dataframe
print("Merging "+str(len(tablelist))+" files into one table")
merged_counts = reduce(lambda left, right: pd.merge(left, right, on = ["Guide"], how = "outer"), tablelist)

# Replace NAs by 0
print("Replacing NAs by 0")
merged_counts = merged_counts.fillna(0)

# Alphabetically reorder dataframe
print("Alphabetically reordering sample columns")
merged_counts = merged_counts.reindex(sorted(merged_counts.columns), axis=1)

# Move Guide column to the first position
print("Getting Guide as the first column")
first_column = merged_counts.pop("Guide")
merged_counts.insert(0, "Guide", first_column)

# Merge counts guide dataframe with controls library to filter only those guides present in the library
print("Filtering by Control Library")
library_file = "/home/jorge/Documents/CRISPRi/library_design/library_controls_final.csv"
library_controls = pd.read_csv(library_file, sep = ";")

# Create a sheet without the gene name for input for downstream code

merged_counts_gene = pd.merge(merged_counts, library_controls, how = "inner", on = ["Guide"])
merged_counts = merged_counts_gene.drop(labels = "name", axis = 1)

# Move gene name to first position

first_column = merged_counts_gene.pop("name")
merged_counts_gene.insert(0, "name", first_column)

# Reformat column names to remove _S*** code from the sequencing demultiplexing

newcols = []
for col in merged_counts.columns:
    if col != "Guide":
        newcol = col.split("_")
        newcol = "_".join(newcol[:-1])
        newcols.append(newcol)
    else:
        newcol = col
        newcols.append(newcol)
    
merged_counts.columns = newcols

newcols = []
for col in merged_counts_gene.columns:
    if col != "Guide" and col != "name":
        newcol = col.split("_")
        newcol = "_".join(newcol[:-1])
        newcols.append(newcol)
    else:
        newcol = col
        newcols.append(newcol)
    
merged_counts_gene.columns = newcols

# Parse second sheet dataframe and get coverage by sample
# By summing the number of reads per sample (column) and multiplying by the number of guides per gene
# Then dividing by the size of the library, we got the approximate coverage per gene by each sample

print("Getting coverage per gene for each sample")
# Initialise coverage row with 2 columns to fix the dimensions to the table
cov_row = ["", "Sample coverage per gene"]
# Declare length of guides library to calculate coverage
lenlibrary = len(merged_counts_gene.index)
for column in merged_counts_gene:
    if column != "Guide" and column != "name":
        mediancountguides = statistics.median(merged_counts_gene[column].values)
        covsample = (mediancountguides*5)
        #print(covsample)
        cov_row.append(covsample)

# Include coverage row in the second sheet
merged_counts_gene.loc[len(merged_counts_gene)] = cov_row
"""
merged_counts_gene = merged_counts_gene.append(pd.DataFrame(cov_row,
                                                            columns = list(merged_counts_gene.columns),
                                                            ignore_index = True))
"""
# Save in an excel df with two identical sheets, the second one containing the gene names

print("Saving xlsx tables")
outname = counts_directory+"merged_counts.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
merged_counts.to_excel(writer, sheet_name='merged counts', index=False)
merged_counts_gene.to_excel(writer, sheet_name='merged counts with gene names', index=False)
writer.save()

print("Output saved to "+counts_directory+outname)
