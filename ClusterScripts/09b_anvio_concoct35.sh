#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_bin_concoct_35_stderr_2020.txt
#SBATCH -o anvio_bin_concoct_35_stderr_2020.txt
#SBATCH -J grass_anvio_bin_concoct
#SBATCH --mem 248000
#SBATCH -t 168:00:00 #time in hours:min:sec
#SBATCH --partition=production



module load anvio/6.1
module load prodigal
module load concoct/1.1.0

source activate concoct-1.1.0
conda activate --stack anvio-6.1

anvi-cluster-contigs -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -C CONCOCT35 --driver concoct -T 48 --clusters 35 --just-do-it
