#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_bin_summarize.txt
#SBATCH -o anvio_bin_summarize.txt
#SBATCH -J grass_anvio_bins
#SBATCH --mem 96G
#SBATCH -t 48:00:00 #time in hours:min:sec
#SBATCH --partition=production


module load anvio/6.2

source activate anvio-6.2

anvi-summarize -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -o GL_meta_merged_kaiju_2500bp/sample_summary_DASTOOL -C DASTOOL

anvi-summarize -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -o GL_meta_merged_kaiju_2500bp/sample_summary_CONCOCT35 -C CONCOCT35

anvi-summarize -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -o GL_meta_merged_kaiju_2500bp/sample_summary_METABAT -C METABAT

anvi-summarize -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -o GL_meta_merged_kaiju_2500bp/sample_summary_BINSANITY -C BINSANITY

anvi-summarize -p GL_meta_merged_kaiju_2500bp/PROFILE.db -c GL.db -o GL_meta_merged_kaiju_2500bp/sample_summary_CONCOCT -C CONCOCT


