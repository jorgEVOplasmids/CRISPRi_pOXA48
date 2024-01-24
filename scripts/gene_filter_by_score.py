#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:34:18 2023

@author: Jorge Sastre Domínguez

This script takes as input the tables summarized by sample and timepoint with log2 scores separated in different
columns. Then it keeps genes with scores above and below +-0.5.

"""

import pandas as pd

results_table = pd.read_excel('/home/jorge/Documents/CRISPRi/results/table_erta_vs_no_erta.xlsx')

# Filter genes with scores below 0.5 both in ERTA and no ERTA conditions

filt_genes = results_table[(results_table["median_log2FC"] > 0.5) | (results_table["median_log2FC"] < -0.5) | (results_table["erta_median_log2FC"] > 0.5) | (results_table["erta_median_log2FC"] < -0.5)]

# Add category column based on the gene score in presence or absence of the antibiotic first separately (1D) then together (2D)

category_col_1D_no_erta = []
category_col_1D_erta = []
category_col_2D = []

for index, row in filt_genes.iterrows():
    score_no_erta = row["median_log2FC"]
    score_erta = row["erta_median_log2FC"]
    if score_no_erta > 0.5:
        category_col_1D_no_erta.append("Costly in absence of ERTA")
    elif score_no_erta < -0.5:
        category_col_1D_no_erta.append("Beneficial in absence of ERTA")
    else:
        category_col_1D_no_erta.append("")
    if score_erta > 0.5:
        category_col_1D_erta.append("Costly in presence of ERTA")
        
    elif score_erta < -0.5:
        category_col_1D_erta.append("Beneficial in presence of ERTA")
    else:
        category_col_1D_erta.append("")

filt_genes["Category LB+DAPG"] = category_col_1D_no_erta
filt_genes["Category LB+DAPG+ERTA"] = category_col_1D_erta

for index, row in filt_genes.iterrows():
    category_no_erta = row["Category LB+DAPG"]
    category_erta = row["Category LB+DAPG+ERTA"]
    if category_no_erta == "Beneficial in absence of ERTA" and category_erta == "Beneficial in presence of ERTA":
        category_col_2D.append("Beneficial for pOXA-48")
    elif category_no_erta == "Costly in absence of ERTA" and category_erta == "Beneficial in presence of ERTA":
        category_col_2D.append("Beneficial under Ab pressure")
    elif category_no_erta == "Costly in absence of ERTA" and category_erta == "Costly in presence of ERTA":
        category_col_2D.append("Costly")
    elif category_no_erta == "Beneficial in absence of ERTA" and category_erta == "Costly in presence of ERTA":
        category_col_2D.append("Costly under Ab pressure")
    elif category_no_erta == "" and category_erta != "":
        category_col_2D.append(category_erta)
    elif category_no_erta != "" and category_erta == "":
        category_col_2D.append(category_no_erta)

filt_genes["Category ERTA vs NO ERTA"] = category_col_2D
filt_genes.drop("gene", axis = 1, inplace = True)
filt_genes.drop("treatment", axis = 1, inplace = True)

outname = "/home/jorge/Documents/CRISPRi/results/gene_scores_with_categories.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
filt_genes.to_excel(writer, sheet_name='scores with categories', index=False)
writer.save()
