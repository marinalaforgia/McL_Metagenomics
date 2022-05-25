#!/bin/bash
#SBATCH -p production
#SBATCH --mem 150G
#SBATCH -n 24
#SBATCH --out Salmon_featureCounts.log
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -J gl_salmon
#SBATCH -t 72:00:00 #time in hours:min:sec

COASSEM=$(realpath megahit_coassembly/GL_coassembly_fixed.fa)
GENE=$(realpath GL_gene_calls.fa)
CPU=24
SAMPFILE=File_names.tsv
COINDEX=GL_salmon_co.idx
INDEX=GL_salmon_gene.idx
OUTDIR=salmon_co
OUTDIR2=salmon_gene

mkdir $OUTDIR
mkdir $OUTDIR2

module load salmon/1.3.0

if [ ! -f $COINDEX ]; then
	salmon index -t $COASSEM -i $COINDEX -k 31

fi

if [ ! -f $INDEX ]; then
	salmon index -t $GENE -i $INDEX -k 31

fi



for prefix in `ls *.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_R1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_R2_filtered.fastq.gz )

	salmon quant -i $INDEX -l IU -1 ${read1} -2 ${read2} --validateMappings -o $OUTDIR2/${prefix}.quant -p $CPU

done


for prefix in `ls *.gz | cut -f1,2 -d'_' | sort -u`;
do
	echo $prefix
	read1=( ${prefix}*_R1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( ${prefix}*_R2_filtered.fastq.gz )

	salmon quant -i $COINDEX -l IU -1 ${read1} -2 ${read2} --validateMappings -o $OUTDIR/${prefix}.quant -p $CPU
	
done



