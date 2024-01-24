for reads in /home/jorge/Documents/CRISPRi/sequences_illumina/AliciaCNB/*.fastq.gz
do
	fastqc $reads -t 7 -o /home/jorge/Documents/CRISPRi/sequences_illumina/fastqc_out
done

multiqc /home/jorge/Documents/CRISPRi/sequences_illumina/fastqc_out -o /home/jorge/Documents/CRISPRi/sequences_illumina/multiqc_out


