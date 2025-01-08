# CRISPRi screening of pOXA-48

## Phenotypic analysis of the strains under study

### Phylogenetic tree

To represent the sequence types considered in our study, we measured the phylogenetic distances between the strains (default [mashtree](https://github.com/lskatz/mashtree) command), and represented the resulting tree using [iTOL](https://itol.embl.de/)), which is shown in the main manuscript (**Figure 1C**).

### pOXA-48 fitness effects

To analyze the fitness effects caused by pOXA-48 in its bacterial hosts, we measured the bacterial growth of each strain both carrying and not carrying pOXA-48. We analyzed these data in R (see SCRIPT). The resulting plot is shown in **Figure 1D**. Furthermore, in following sections of the paper, we analyze whether the deletion of the resistance gene from the plasmid alleviates the plasmid costs. These analyses can be found in the script [growth_curves_delta_blaOXA48.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/growth_curves_delta_blaOXA48.R).

## Analysis of CRISPRi screen

#### Quality Control of raw reads

First, we performed a QC of the raw Illumina data sequences (55 bp amplicon sequencing) using FastQC and MultiQC. After confirming high quality per base score in the guides region (mean Phred score > 30), we continued with the downstream analyses (see [qc_crispr.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/qc_crispr.sh)).

#### Guide detection, count and reformatting

For getting each individual guide count in each sample, we used Bash, and parsed each fastq file. Knowing the sequence of the promoter preceeding the guides in each read, we used 15 bp to match the reads and extract the following 20 bp (gRNA). With bash commands we got for each sample the number of times each individual guide appeared. We sorted the list of guides, counted each unique occurrence, and reformatted the output to save each file count in an individual csv (see [counts_script.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/counts_script.sh)).

We then took all the individual csv counts files and filtered the guides with less counts than 5 fot timepoints 0 and with less counts than 2 for the rest of the timepoints. We merged all the samples in the same table, and filtered only those guides present in the plasmid library. We got as output an excel called merged_counts.xlsx, with two sheets, one with the gene names and another without them, for following code execution (see [merge_counts.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/merge_counts.py)). In the second sheet, we also include the sequencing coverage information for each demultiplexed sample, approximated as the sequencing coverage per gene (which is desired to be > 300x) calculating it as the median of the counts of guides (1 gRNA per read; avoiding outliers from guides counted in excess) multiplied by 5 (as there are 5 guides per gene).

We then took the counts of the guides table by replicate and merged the counts by sample (summed the three individual replicates by timepoint per each strain (see [CRISPR_merge_replicates.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/CRISPR_merge_replicates.py)) after checking the reproducibility of the replicates in each condition (LB+DAPG or LB+DAPG+ERTA). After this step we obtained the table which will be the input for the downstream code.

#### Gene score calculation

For calculating the gene scores (median log2 fold change by gene), again, based on the code of F. Rousset ([https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb](https://gitlab.pasteur.fr/dbikard/ecocg/-/blob/master/Notebook.ipynb)) we wrote the [CRISPR_pool_log2FC_scores.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/CRISPR_pool_log2fc_scores.py) script, which takes the table with the replicates merged by sample, and calculates the log2FC for each gene, then centering the values by removing the median of control guides. The output table contains two sheets, being the first one the table with log2FC values and the second one the table with the metadata.

#### Preeliminary plotting

Finally, to represent the resulting gene scores in an easy way, we decided to add the annotation of the plasmid genes and intergenic regions (see [pOXA_intergenic_reannotation.py](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/pOXA_intergenic_reannotation.py)). Then, we represented the log2FC for each gene in each timepoint of each treatment for each strain, colouring by annotation function. These first graphics give an idea of when is it beneficial to silence a gene (positive median gene score) or deleterium when silencing it (negative median gene score), with the [plots_log2FC_by_gene.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/plots_log2FC_by_gene.R) script. We did the same representing in the Y axis the ERTA treatment and in the X axis the NO ERTA treatment.

#### Final plot of CRISPRi scores

To represent the complete information including all the strains after processing the data, we decided to represent the gene scores in a heatmap for each experimental condition (**Figure 2A**). These were plotted using R (see [heatmap_CRISPRi.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/heatmap_CRISPRi.R)).

## Analysis of response conservation

To analyze whether the response to the CRISPRi screens in the different conditions was conserved among the diverse strains, we correlated the gene scores for each pairwise strain combination per condition using R corrplot package (see [similarity_matrices_and_SDC_analysis.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/similarity_matrices_and_SDC_analysis.R)). These correlation matrices are shown in **Figure 2B**. Additionally, this script contains the code to perform the similarity analysis using the Sørensen-Dice Coefficient (SDC), which are shown in **Supplementary Figure 5**.

## Permutation test

For checking statistically significant differences from the score of each gene from the distribution of gene scores, we performed a permutation test for all the data merged (14 strains & 106 genes). Taking each gene score, we calculated the difference between this value and the mean of the rest of the gene scores. Then, for getting the probability of obtaining the score of each gene, we calculated for 100000 permutations  the difference between random gene scores and the mean of the rest of the genes (see [CRISPRi_stats.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/CRISPR_stats.R)). Lastly, the proportion of times that this statistic was more extreme than the score of each gene indicated the probability of getting such a score (i.e. 14 measurements per gene, calc its mean, compare vs rest of 100000 random means calculated, minimum p-value = 0.00001). We selected 100000 permutations to get more precise probabilities, as the p values of some genes seemed marginally significant when performing the permutation tests with 10000 repeats. Finally, we adjusted the p-values by FDR. The statistical analyses of the CRISPRi results through the application of these permutation tests are shown in the **Figure 3A**.

## Correlation between CRISPRi scores and pOXA-48 costs

To study the relationship between the blaOXA-48 gene score and the plasmid cost, we calculated the correlation between both variables (see [correlations_costs_CRISPRi_scores.R](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/correlations_costs_CRISPRi_scores.R)), and represented it as shown in **Figure 5A**.
