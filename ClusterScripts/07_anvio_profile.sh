#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_stderr.txt
#SBATCH -o anvio_stderr.txt
#SBATCH -J grass_anvio_profile
#SBATCH --mem 490000
#SBATCH -t 169:00:00 #time in hours:min:sec
#SBATCH --partition=production


module load anvio/6.1
module load prodigal

source activate anvio-6.1




for sample in $(cat File_names_v2020A.txt); 
do anvi-profile -i megahit_coassembly_mapping/$sample'.bam' -c GL.db --num-threads 24 --min-contig-length 2500; 
done



