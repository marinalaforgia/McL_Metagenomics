#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_bin_binsanity_stderr.txt
#SBATCH -o /anvio_bin_binsanity_stderr.txt
#SBATCH -J grass_anvio_bin_binsanity
#SBATCH --mail-type=END # notifications for job done & fail
#SBATCH --mail-user=clettinger@ucdavis.edu # send-to address
#SBATCH --mem 248000
#SBATCH -t 168:00:00 #time in hours:min:sec
#SBATCH --partition=production



module load anvio/6.1
module load binsanity/0.3.3 

source activate anvio-6.1

anvi-cluster-contigs -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -C BINSANITY --driver binsanity -T 48 --just-do-it
