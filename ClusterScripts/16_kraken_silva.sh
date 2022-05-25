#!/bin/bash
#SBATCH -p production
#SBATCH --mem 96G
#SBATCH -n 16
#SBATCH --out kraken_silva.log
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -J gl_kraken_silva
#SBATCH -t 72:00:00 #time in hours:min:sec


module load kraken2/2.1.1
module load bracken/2.5


CPU=16

export KRAKEN2_DB_PATH="/share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/"

bracken-build -d silva -t $CPU -k 35 -l 150

for prefix in `ls *.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_R1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_R2_filtered.fastq.gz )

	kraken2 --db silva --threads $CPU --report kraken2_report_paired_silva_${prefix}.tsv --gzip-compressed --paired ${read1} ${read2} > ${prefix}.silva.log

	bracken -d silva -i kraken2_report_paired_silva_${prefix}.tsv -o ${prefix}.silva.S.bracken -r 150 -l 'S' -t 10
	bracken -d silva -i kraken2_report_paired_silva_${prefix}.tsv -o ${prefix}.silva.F.bracken -r 150 -l 'F' -t 10


done

