# CRISPRi screening of pOXA-48

## Analysis of CRISPRi screen

#### Quality Control of raw reads

First, we performed a QC of the raw Illumina data sequences (55 bp amplicon sequencing) using FastQC and MultiQC. After confirming high quality per base score in the guides region (mean Phred score > 30), we continued with the downstream analyses (see [qc_crispr.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/qc_crispr.sh)).

#### Guide detection, count and reformatting

For getting each individual guide count in each sample, we used Bash, and parsed each fastq file. Knowing the sequence of the promoter preceeding the guides in each read, we used 15 bp to match the reads and extract the following 20 bp (gRNA). With bash commands we got for each sample the number of times each individual guide appeared. We sorted the list of guides, counted each unique occurrence, and reformatted the output to save each file count in an individual csv (see [counts_script.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/counts_script.sh)).

We then took all the individual csv counts files and filtered the guides with less counts than 5 fot timepoints 0 and with less counts than 2 for the rest of the timepoints. We merged all the samples in the same table, and filtered only those guides present in the plasmid library. We got as output an excel called merged_counts.xlsx, with two sheets, one with the gene names and another without them, for following code execution (see [merge_counts.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/merge_counts.py)). In the second sheet, we also include the sequencing coverage information for each demultiplexed sample, approximated as the sequencing coverage per gene (which is desired to be > 300x) calculating it as the median of the counts of guides (1 gRNA per read; avoiding outliers from guides counted in excess) multiplied by 5 (as there are 5 guides per gene).

We then took the counts of the guides table by replicate and merged the counts by sample (summed the three individual replicates by timepoint per each strain (see [CRISPR_merge_replicates.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/CRISPR_merge_replicates.py)) after checking the reproducibility of the replicates in each condition (LB+DAPG or LB+DAPG+ERTA). After this step we obtained the table which will be the input for the downstream code.

#### Gene score calculation

For calculating the gene scores (median log2 fold change by gene), again, based on the code of F. Rousset ([https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb](https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb)) we wrote the [CRISPR_pool_log2FC_scores.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/CRISPR_pool_log2fc_scores.py) script, which takes the table with the replicates merged by sample, and calculates the log2FC for each gene, then centering the values by removing the median of control guides. The output table contains two sheets, being the first one the table with log2FC values and the second one the table with the metadata.

#### Preeliminary plotting and filtering by score

Finally, to represent the resulting gene scores in an easy way, we decided to add the annotation of the plasmid genes and intergenic regions (see [pOXA_intergenic_reannotation.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/pOXA_intergenic_reannotation.py)). Then, we represented the log2FC for each gene in each timepoint of each treatment for each strain, colouring by annotation function. These first graphics give an idea of when is it beneficial to silence a gene (positive median gene score) or deleterium when silencing it (negative median gene score), with the [plots_log2FC_by_gene.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/plots_log2FC_by_gene.R) script. We did the same representing in the Y axis the ERTA treatment and in the X axis the NO ERTA treatment.

With these data, we summarized the genes in tables depending on the values of the log2FC, setting as threshold -0.5 and +0.5, and added category columns for each condition and also combining both conditions (ERTA vs no ERTA). We did this with the [gene_filter_by_score.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/gene_filter_by_score.py) script.
