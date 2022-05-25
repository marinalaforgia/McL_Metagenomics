#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH -J GL_prot_annot
#SBATCH --mem=200G #memory
#SBATCH -p intel,batch,stajichlab
#SBATCH -D /rhome/cassande/bigdata/eisenlab/grassland_MAGs/
#SBATCH -o logs/13b_annot.log
#SBATCH -e logs/13b_annot.log


module unload miniconda2
module unload anaconda3
module load miniconda3
module load interproscan/5.52-86.0
module load diamond

source activate anvio-7

anvi-get-sequences-for-gene-calls -c GL.db --get-aa-sequences -o GL_amino_acid_sequences.fa

source deactivate


module load eggnog-mapper/1.0.3

emapper.py -i GL_amino_acid_sequences.fa --output GL_prot -m diamond --cpu 24


source activate anvio-7

anvi-script-run-eggnog-mapper -c GL.db --annotation GL_prot.emapper.fix.annotations --use-version 1.0.3

anvi-delete-hmms -c GL.db --just-do-it
anvi-run-hmms -c GL.db --num-threads 24 --just-do-it

anvi-run-pfams -c GL.db --num-threads 24
anvi-run-kegg-kofams -c GL.db --num-threads 24

anvi-run-ncbi-cogs -c GL.db --num-threads 24



module load interproscan/5.52-86.0
interproscan.sh -i GL_amino_acid_sequences.fa -f tsv -o GL_interpro-output.tsv


anvi-import-functions -c GL.db -i GL_interpro-output.tsv -p interproscan


anvi-get-sequences-for-gene-calls -c GL.db -o GL_AA_v2.gff --export-gff3

