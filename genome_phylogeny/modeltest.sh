#!/bin/bash
#PBS -k o 
#PBS -l nodes=1:ppn=8,vmem=32gb,walltime=36:00:00
#PBS -M cholden23@gmail.com
#PBS -m abe 
#PBS -N modeltest
#PBS -j oe 


#####################################################
#####################################################

cd /N/dc2/scratch/chrichan/wolbachia/phylogeny

modeltest-ng -d nt \
 -i wolbachia_polymorphisms.fasta \
 -o modeltest_output.txt \
 -p 8 \
 -h uigf \
 -f ef 
 