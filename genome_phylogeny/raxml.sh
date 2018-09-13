#!/bin/bash
#PBS -k o 
#PBS -l nodes=1:ppn=8,vmem=32gb,walltime=36:00:00
#PBS -M cholden23@gmail.com
#PBS -m abe 
#PBS -N raxml
#PBS -j oe 


#####################################################
#####################################################

module load raxml/gnu/8.2.11

cd /N/dc2/scratch/chrichan/wolbachia/phylogeny

#Rapid run
raxmlHPC-PTHREADS -f a \
 -m GTRGAMMAIX \
 -p 39812 \
 -x 11253 \
 -# 200 \
 -s wolbachia_polymorphisms.fasta \
 -n rapid_bootstrap \
 -T 8

#Do an ML search for best-scoring tree
raxmlHPC-PTHREADS -s wolbachia_polymorphisms.fasta \
 -n mlsearch \
 -m GTRGAMMAIX \
 -p 54376 \
 -T 8 \
 -N 100 

#Next, do a bootstrap search
raxmlHPC-PTHREADS -s wolbachia_polymorphisms.fasta \
 -n bootstrapsearch \
 -m GTRGAMMAIX \
 -b 32198 \
 -p 87632 \
 -T 8 \
 -N 200 

#Now assign bootstrap support values to nodes

raxmlHPC-PTHREADS  -m GTRGAMMAIX \
 -p 41294 \
 -f b \
 -t RAxML_bestTree.mlsearch \
 -z RAxML_bootstrap.bootstrapsearch \
 -n mlwithbootstrap \
 -T 8




#raxmlHPC-SSE3 -s wolbachia_polymorphisms.fasta -m GTRGAMMAIX -n EXEC_NAME -p PARSIMONY_SEED
#raxml-ng --msa wolbachia_polymorphisms.fasta --model TPM1uf+I+G4

#raxmlHPC-SSE3 -s wolbachia_polymorphisms.fasta -m GTRGAMMAIX -n EXEC_NAME -p PARSIMONY_SEED
#raxml-ng --msa wolbachia_polymorphisms.fasta --model TVM+I+G4
