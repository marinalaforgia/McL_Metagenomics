#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e bam_stderr.txt
#SBATCH -o bam_stderr.txt
#SBATCH -J grass_coverage
#SBATCH --mem 98000
#SBATCH -t 144:00:00 #time in hours:min:sec
#SBATCH --partition=production

module load bowtie2
module load samtools
module load anvio

source activate anvio5


for file in $(cat File_names.txt); 
do bowtie2 --threads 16 -x megahit_coassembly_mapping/GL_contigs -1 $file'_R1_filtered.fastq' -2 $file'_R2_filtered.fastq' -S megahit_coassembly_mapping/$file'.sam'; 
done

for file in $(cat File_names.txt); 
do samtools view -F 4 -bS megahit_coassembly_mapping/$file'.sam' > megahit_coassembly_mapping/$file'-RAW.bam'; 
done

for file in $(cat File_names.txt); 
do anvi-init-bam megahit_coassembly_mapping/$file'-RAW.bam' -o megahit_coassembly_mapping/$file'.bam'; 
done

