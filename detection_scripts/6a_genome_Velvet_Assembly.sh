#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -N assembly_velvet
#PBS -j oe

# This script will use the assembler Velvet to assemble the genes that have matches with the Wolbachia genes
# For each assembly multiple kmer sizes are tested
# The fasta files that were produced from the make_wol_fasta.sh script is the input
# The output will be stored in the ~/assembly_velvet directory


module load velvet

# ftsZ genes
# Make the accession numbers of positive samples a variable
CUR_ACC=~/downloads/positive_samples_subset.txt

while read CUR_ACC; do
  echo ${CUR_ACC}
  for kmer in 21 31 41 51
  do
    velveth ~/assembly_velvet/${CUR_ACC}_ftsZ_${kmer} ${kmer} ~/positive_fasta_files/ftsZ/${CUR_ACC}_ftsZ.fasta
    velvetg ~/assembly_velvet/${CUR_ACC}_ftsZ_${kmer} -cov_cutoff auto
  done
done <${CUR_ACC}
echo 'done with ftsZ ***********************************************************************************************'

# groE genes
CUR_ACC=~/downloads/positive_samples_subset.txt
while read CUR_ACC; do
  echo ${CUR_ACC}
  for kmer in 21 31 41 51
  do
    velveth ~/assembly_velvet/${CUR_ACC}_groE_${kmer} ${kmer} ~/positive_fasta_files/groE/${CUR_ACC}_groE.fasta
    velvetg ~/assembly_velvet/${CUR_ACC}_groE_${kmer} -cov_cutoff auto
  done
done <${CUR_ACC}
echo 'done with groE ***********************************************************************************************'

# wsp genes
CUR_ACC=~/downloads/positive_samples_subset.txt
while read CUR_ACC; do
  echo ${CUR_ACC}
  for kmer in 21 31 41 51
  do  
    velveth ~/assembly_velvet/${CUR_ACC}_wsp_${kmer} ${kmer} ~/positive_fasta_files/wsp/${CUR_ACC}_wsp.fasta
    velvetg ~/assembly_velvet/${CUR_ACC}_wsp_${kmer} -cov_cutoff auto
  done
done <${CUR_ACC}
echo 'done with wsp ***********************************************************************************************'

mv ~/assembly_velvet/*wsp* ~/assembly_velvet/wsp/
mv ~/assembly_velvet/*groE* ~/assembly_velvet/groE/
mv ~/assembly_velvet/*ftsZ* ~/assembly_velvet/ftsZ/