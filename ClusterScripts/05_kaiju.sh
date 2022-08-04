#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/QC_fastqs/megahit_coassembly/
#SBATCH -e cassie_kaiju_stderr.txt
#SBATCH -J grassland_co_kaiju
#SBATCH --mem-per-cpu=4G #memory per node in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec
#SBATCH --partition=production

/share/eisenlab/gjospin/software/kaiju/src/kaiju -z 24 -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -f /share/eisenlab/casett/database/kaiju_nr_euk/kaiju_db_nr_euk.fmi -i final.contigs.fa -o kaiju.out.updated

/share/eisenlab/gjospin/software/kaiju/src/kaijuReport -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i kaiju.out.updated -r genus -p -o kaiju.out.summary.updated

/share/eisenlab/gjospin/software/kaiju/src/addTaxonNames -t /share/eisenlab/casett/database/kaiju_nr_euk/nodes.dmp -n /share/eisenlab/casett/database/kaiju_nr_euk/names.dmp -i kaiju.out.updated -o kaiju.names.out.updated.samelevels -r phylum,order,class,family,genus,species
