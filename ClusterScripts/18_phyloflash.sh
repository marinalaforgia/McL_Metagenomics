#!/bin/bash
#SBATCH -p production
#SBATCH --mem 148G
#SBATCH -n 24
#SBATCH --out phyloflash.log
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -J gl_phyloflash
#SBATCH -t 72:00:00 #time in hours:min:sec

module load phyloflash

CPU=24

#phyloFlash_makedb.pl --remote

for prefix in `ls *.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_R1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_R2_filtered.fastq.gz )

	phyloFlash.pl -lib ${prefix} -almosteverything -CPUs $CPU -readlength 150 -read1 ${read1} -read2 ${read2}

done


