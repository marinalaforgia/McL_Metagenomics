#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/
#SBATCH -e GL_kaiju_stderr.txt
#SBATCH -J grassland_kaiju_import
#SBATCH --mem-per-cpu=4G #memory per node in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec
#SBATCH --partition=production


module load anvio/6.2
source activate anvio-6.2

anvi-get-sequences-for-gene-calls -c GL.db -o GL_gene_calls.fa

/share/eisenlab/gjospin/software/kaiju/src/kaiju -z 24 -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -f /share/eisenlab/casett/database/kaiju_nr_euk/kaiju_db_nr_euk.fmi -i GL_gene_calls.fa -o GL_gene_calls.kaiju.out.updated -v

/share/eisenlab/gjospin/software/kaiju/src/kaijuReport -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i GL_gene_calls.kaiju.out.updated -r genus -p -o GL_gene_calls.summary.updated

/share/eisenlab/gjospin/software/kaiju/src/addTaxonNames -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i GL_gene_calls.kaiju.out.updated -o GL_gene_calls.kaiju.names.out.updated.samelevels_v2 -r superkingdom,phylum,class,order,family,genus,species

anvi-import-taxonomy-for-genes -i GL_gene_calls.kaiju.names.out.updated.samelevels_v2 -c GL.db -p kaiju --just-do-it
