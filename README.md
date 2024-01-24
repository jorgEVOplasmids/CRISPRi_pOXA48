# 

## Analysis of CRISPRi screen

#### Quality Control of raw reads

First, we performed a QC of the raw Illumina data sequences (55 bp amplicon sequencing) using FastQC and MultiQC. After confirming high quality per base score in the guides region (mean Phred score > 30), we continued with the downstream analyses (see [qc_crispr.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/qc_crispr.sh)).

#### Guide detection, count and 

For getting each individual guide count in each sample, we used Bash, and parsed each fastq file. Knowing the sequence of the promoter preceeding the guides in each read, we used 15 bp to match the reads and extract the following 20 bp (gRNA). With bash commands we got for each sample the number of times each individual guide appeared. We sorted the list of guides, counted each unique occurrence, and reformatted the output to save each file count in an individual csv (see [counts_script.sh](https://github.com/jorgEVOplasmids/CRISPRi_pOXA48/blob/main/scripts/counts_script.sh)).


