#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_merge_stderr.txt
#SBATCH -o anvio_merge_stderr.txt
#SBATCH -J grass_anvio_merge
#SBATCH --mem 490000
#SBATCH -t 168:00:00 #time in hours:min:sec
#SBATCH --partition=production



module load anvio/6.1
module load prodigal

source activate anvio-6.1



anvi-merge megahit_coassembly_mapping/*ANVIO_PROFILE/PROFILE.db -o GL_meta_merged_kaiju_2500bp/ -c GL.db --enforce-hierarchical-clustering

