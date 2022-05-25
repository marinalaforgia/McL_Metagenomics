#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_bin_metabat_stderr.txt
#SBATCH -o anvio_bin_metabat_stderr.txt
#SBATCH -J grass_anvio_bin_metabat
#SBATCH --mem 248000
#SBATCH -t 168:00:00 #time in hours:min:sec
#SBATCH --partition=production



module load anvio/6.1
module load prodigal
module load metabat/2.12.1

source activate anvio-6.1

anvi-cluster-contigs -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -C METABAT --driver metabat2 -T 48 --just-do-it
