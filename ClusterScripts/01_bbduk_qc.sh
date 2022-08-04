#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -D /share/eisenlab/casett/grassland/raw_fastqs/combined_fastq/
#SBATCH -e bbduk_stderr.txt
#SBATCH -o bbduk_stderr.txt
#SBATCH -J grassland_qc
#SBATCH --partition=production
#SBATCH --mem-per-cpu=2G #memory per node in Gb
#SBATCH -t 48:00:00 #time in hours:min:sec

module load bbmap/37.68

for file in $(cat File_names.txt); 
do bbduk.sh in1=$file'_R1_combined.fastq' in2=$file'_R2_combined.fastq' out1=$file'_R1_filtered.fastq' out2=$file'_R2_filtered.fastq' qtrim=rl trimq=10 maq=10; 
done

