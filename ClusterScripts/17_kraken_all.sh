#!/bin/bash
#SBATCH -p production
#SBATCH --mem 248G
#SBATCH -n 32
#SBATCH --out kraken_all.log
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -J gl_kraken_all
#SBATCH -t 72:00:00 #time in hours:min:sec

module load kraken2/2.1.1
module load bracken/2.5

CPU=32

export KRAKEN2_DB_PATH="/share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/"

#kraken2-build --standard --db stan --threads 32

#wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz

#tar -zxvf k2_standard_20210517.tar.gz -C stan


#bracken-build -d stan -t $CPU -k 35 -l 150

for prefix in `ls *.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_R1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_R2_filtered.fastq.gz )

	kraken2 --db stan --threads $CPU --report kraken2_report_paired_${prefix}.tsv --gzip-compressed --paired ${read1} ${read2} > ${prefix}.log

	bracken -d stan -i kraken2_report_paired_${prefix}.tsv -o ${prefix}.S.bracken -r 150 -l 'S' -t 10
	bracken -d stan -i kraken2_report_paired_${prefix}.tsv -o ${prefix}.F.bracken -r 150 -l 'F' -t 10


done

