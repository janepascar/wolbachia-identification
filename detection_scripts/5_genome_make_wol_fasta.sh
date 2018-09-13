#!/bin/bash

#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -N make_wol_fasta
#PBS -j oe

# This script is used to convert the sam files produced from Magic-BLAST into fasta files
# The R script will determine which gene the read matches based on the dictionary.txt file
# A fasta file will be created for each gene for each accession number

# Set the subset of positive accession numbers as a variable
ACC=~/downloads/positive_samples_subset.txt

module load RMod

# Start a loop that runs through all of the accession numbers 
while read ACC; do
  echo ${ACC}
  Rscript --vanilla ~/5_make_wol_fasta.R ~/dictionary.txt ~/Magic_BLAST_results/sam_fmt/${ACC}.sam ~/positive_fasta_files/${ACC}
done <${ACC}

# Organize the fasta files by gene
cd ~/positive_fasta_files
mv *wsp.fasta ~/positive_fasta_files/wsp
mv *ftsZ.fasta ~/positive_fasta_files/ftsZ
mv *groE.fasta ~/positive_fasta_files/groE
