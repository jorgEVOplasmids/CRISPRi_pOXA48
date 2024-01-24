#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 16:52:41 2023

@author: Jorge Sastre Domínguez

This script takes the annotation in gbk format of pOXA-48 and reannotates intergenic regions depending on the
surrounding genes and direction of the intergenic region.

"""

from Bio import SeqIO
import pandas as pd

# GBK file from which we are going to reannotate the intergenic regions to the surrounding genes
gbk_dir = "/home/jorge/Documents/CRISPRi/library_design/annotations/pOXA48K8_JAVI_intergenic_twoorientations.gbk"
# Declare an empty list which willl contain the original gene identifiers to then have a table with both nomenclatures
annot_original = []
# And an empty list for the new names
annot_new = []
# Parse gbk file for each feature
with open(gbk_dir) as gbkfile:
    for record in SeqIO.parse(gbkfile, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                annot_original.append(feature.qualifiers["label"][0]) # Include all the names into the original annotation list
                                   
#print(annot_original)
# From this list, we can rename the intergenic regions, as it is an ordered list of the plasmid genes
for i in range(0,len(annot_original)-2): # Because last 2 elements are intergenic regions and will give index range errors
    #print(annot_new)
    gene_downstream = ""
    gene_upstream = ""
    # Avoid error with last intergenic region in list of genes
    if "intergenic" in annot_original[i] or "intragenic" in annot_original[i]:
        orientation = annot_original[i][-1]
        #print(orientation)
        if "intergenic" in annot_original[i-1] or "intragenic" in annot_original[i-1]:
            gene_downstream = annot_original[i-2]
        else:
            gene_downstream = annot_original[i-1]
        if "intergenic" in annot_original[i+1] or "intragenic" in annot_original[i+1]:
            gene_upstream = annot_original[i+2]
        else:
            gene_upstream = annot_original[i+1]
        annot_new.append(gene_downstream+"/"+gene_upstream+"_"+orientation)
    else:
        annot_new.append(annot_original[i])
        
# Add last intergenic regions manually

annot_new.append(annot_original[-3]+"/"+annot_original[0]+"_A")
annot_new.append(annot_original[-3]+"/"+annot_original[0]+"_B")

# Save in a dictionary and convert into table to export

dic_newannot = {"gene": annot_original, "annot_gene": annot_new}
table_new_annot = pd.DataFrame(dic_newannot)
print(list(table_new_annot["annot_gene"]))
"""
outname = "/home/jorge/Documents/CRISPRi/new_annot_pOXA48.xlsx"
writer = pd.ExcelWriter(outname, engine='xlsxwriter')
table_new_annot.to_excel(writer, sheet_name='complete annotations', index=False)
writer.save()
"""