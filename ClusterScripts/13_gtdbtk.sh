#!/bin/bash -l
#SBATCH --ntasks=16 # Number of cores
#SBATCH --mem=400G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p intel,batch
#SBATCH -o logs/12_gtdbtk.log
#SBATCH -e logs/12_gtdbtk.log
#SBATCH -J grass_gtdbtk
#SBATCH -D /rhome/cassande/bigdata/eisenlab/grassland_MAGs/


module unload miniconda2
module load miniconda3

conda activate gtdbtk-1.5.0

INPUT=bin_fasta

OUTPUT=gtbdk_results
CPU=16
PREFIX=GL

gtdbtk classify_wf --genome_dir $INPUT --out_dir $OUTPUT -x .fa --cpus $CPU --prefix $PREFIX.gtbdk
