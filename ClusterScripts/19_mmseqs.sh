#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH -o logs/mmeseqdb.log
#SBATCH -e logs/mmeseqdb.log
#SBATCH --mem=200G #memory
#SBATCH -J GL_mmseq

source activate mmseqs2

#mkdir tmp

mmseqs databases NR NRref tmp

mmseqs createtaxdb NRref tmp

mmseqs createindex NRref tmp

mmseqs easy-taxonomy GL_amino_acid_sequences_EUK.fa NRref taxonomyResult_nr_GL_euk2 tmp2 --lca-ranks kingdom,phylum,class,order,family,genus --tax-lineage 1
mmseqs easy-taxonomy GL_amino_acid_sequences_UNCL.fa NRref taxonomyResult_nr_GL_uncl2 tmp2 --lca-ranks kingdom,phylum,class,order,family,genus --tax-lineage 1

mmseqs easy-search GL_amino_acid_sequences_EUK.fa NRref functionResult_nr_GL_euk tmp2
mmseqs easy-search GL_amino_acid_sequences_UNCL.fa NRref functionResult_nr_GL_uncl tmp2

