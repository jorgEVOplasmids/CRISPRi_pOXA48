for reads in /home/jorge/Documents/CRISPRi/sequences_illumina/AliciaCNB/*.fastq.gz
do
	samplename=$( basename $reads )
	samplename=${samplename%%_L001_R1_001.fastq.gz}
	#echo $samplename
	echo $samplename" Guide">> /home/jorge/Documents/CRISPRi/guides_count/counts_$samplename.csv
	zcat $reads | grep -o "AGGTATAATACTAGT.\{20\}" | cut -c 16-35 | sort | uniq -c | tr -s " " | sed 's/ //' >> /home/jorge/Documents/CRISPRi/guides_count/counts_$samplename.csv
	echo "Guides counted for " $samplename
done
