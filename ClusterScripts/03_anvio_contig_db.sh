#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e anvio_stderr_2020.txt
#SBATCH -o anvio_stderr_2020.txt
#SBATCH -J grass_anvio_db
#SBATCH --mem 200000
#SBATCH -t 144:00:00 #time in hours:min:sec
#SBATCH --partition=production

module load anvio/6.1
module load prodigal

source activate anvio-6.1

anvi-gen-contigs-database -f megahit_coassembly/GL_coassembly_fixed.fa -o GL.db

anvi-run-hmms -c GL.db --num-threads 24


anvi-run-hmms -T 24 -c GL.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Bacteria_71
anvi-run-hmms -T 24 -c GL.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Archaea_76
anvi-run-hmms -T 24 -c GL.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Ribosomal_RNAs
anvi-run-hmms -T 24 -c GL.db -H /share/eisenlab/gjospin/software/anvio/anvio/data/hmm/Protista_83


anvi-get-sequences-for-gene-calls -c GL.db -o GL_gene_calls.fa




