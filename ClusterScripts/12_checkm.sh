#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -p intel,batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/grassland_MAGs/
#SBATCH -o logs/11_checkm.log
#SBATCH -e logs/11_checkm.log
#SBATCH -J grass_checkM
#SBATCH --mem 96G #memory in Gb
#SBATCH -t 12:00:00 #time in hours:min:sec


module load checkm

BINFOLDER=bin_fasta
OUTPUT=dastool_checkM
CPU=24

checkm lineage_wf -t $CPU -x fa $BINFOLDER $OUTPUT

checkm tree $BINFOLDER -x .fa -t $CPU $OUTPUT/tree

checkm tree_qa $OUTPUT/tree -f $OUTPUT/$OUTPUT.checkm.txt
