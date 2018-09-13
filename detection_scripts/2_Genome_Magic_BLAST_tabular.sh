#!/bin/bash

#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N 2_test_genome
#PBS -j oe

# This script will use Magic-BLAST to map sequences from the SRA to the Wolbachia genes selected
# Three file formats will be created: tabular, csv, and sam

# Change into the directory where the results of Magic-BLAST will be stored
cd ~/Magic_BLAST_results

# Make the list of accession numbers a variable
ACC=~/downloads/accession_numbers.txt

# Make the list of reference Wolbachia genes a variable
REF_GENES=~/Magic_BLAST_results/reference_wol_genes.fasta

# Make a database using the selected reference genes
makeblastdb -in ${REF_GENES} -dbtype nucl

# Run a loop to blast all downloaded data
# Output here will be in tabular format
while read CUR_ACC; do
  echo ${CUR_ACC}
  magicblast -infmt fastq -query ~/downloads/SRA_sequences/${CUR_ACC}.fastq -db ${REF_GENES} -out ~/Magic_BLAST_results/tabular_fmt/${CUR_ACC}.blastresults -outfmt tabular -num_threads 8
  # For downstream analysis files must be in .csv format
  cut -f 1,2,3,7,8,9,10,13,25 ~/Magic_BLAST_results/tabular_fmt/${CUR_ACC}.blastresults >~/Magic_BLAST_results/csv_fmt/${CUR_ACC}.csv
done <${ACC}

echo 'done with Magic-BLAST *********************************************************************************'
