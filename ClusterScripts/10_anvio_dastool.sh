#!/bin/bash -l
#
#SBATCH -n 48 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_bin_dastools_stderr.txt
#SBATCH -o anvio_bin_dastools_stderr.txt
#SBATCH -J grass_anvio_bin_das
#SBATCH --mem 248000
#SBATCH -t 168:00:00 #time in hours:min:sec
#SBATCH --partition=production

export R_LIBS="$R_LIBS:/share/eisenlab/gjospin/R_libs/"

module load anvio/6.2
module load prodigal
module load dastool/1.1.2
module load diamond/0.9.24

source activate dastool-1.1.2
conda activate --stack anvio-6.2

anvi-cluster-contigs -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -C DASTOOL --driver dastool -S CONCOCT,BINSANITY,METABAT,CONCOCT35 -T 48 --search-engine diamond --just-do-it
