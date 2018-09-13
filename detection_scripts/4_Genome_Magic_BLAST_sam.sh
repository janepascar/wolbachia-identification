#!/bin/bash

#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=48:00:00
#PBS -N Magic_BLAST_sam
#PBS -j oe

# Once a list of positive samples based on the presence of wsp genes the samples will be Magic_BLAST again
# This time the references genes will be wsp, FtsZ, and groE genes
# A maximum of 3 samples from each species will be used for phylogenetic analysis

# Set the accession numbers and the reference genes as variables
ACC=~/downloads/positive_samples_subset.txt
REF_GENES=~/Magic_BLAST_results/reference_wol_genes.fasta

# Create a loop to blast all of the downloaded data against the reference genes
# The output will be to sam format
makeblastdb -in ${REF_GENES} -dbtype nucl
while read CUR_ACC; do
  echo ${CUR_ACC}
  magicblast -infmt fastq -query ~/downloads/SRA_sequences/${CUR_ACC}.fastq -db ${REF_GENES} -out ~/Magic_BLAST_results/sam_fmt/${CUR_ACC}.sam -outfmt sam -num_threads 8
done <${ACC}

echo 'done with sam format *********************************************************************************'

# Once blasting is complete the SRA fastq files can be deleted so they are not monopolizing cluster space