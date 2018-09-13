#!/bin/bash
#PBS -k o 
#PBS -l nodes=1:ppn=8,vmem=32gb,walltime=48:00:00
#PBS -M cholden23@gmail.com
#PBS -m abe 
#PBS -N XXSPECIESXX
#PBS -j oe 


#####################################################
#####################################################

SPECIES=XXSPECIESXX
WORKDIR=/N/dc2/scratch/chrichan/wolbachia/${SPECIES}
ACCESSIONS='XXACCESSIONXX'
REF=/N/dc2/scratch/chrichan/wolbachia/assembled_genomes/${SPECIES}.fasta

mkdir -p ${WORKDIR}
cd ${WORKDIR}

###############################
#Download SRA data
mkdir -p raw_data
cd raw_data
for X in ${ACCESSIONS}
do
    fastq-dump --split-files ${X}
done

cat *_1.fastq > all_F.fastq
cat *_2.fastq > all_R.fastq

cd ..

###############################
# Map reads to genome using bwa
bwa index ${REF}

bwa mem -t 8 ${REF} raw_data/all_F.fastq raw_data/all_R.fastq | sambamba view -S -f bam -o ${SPECIES}_unsorted.bam /dev/stdin

sambamba sort -p -t 8 -m 24G -o ${SPECIES}.bam ${SPECIES}_unsorted.bam

###############################
# Generate coverage

#Get length info for each contig in the assembly
seqtk comp ${REF} > seq_info.txt

#Generate list of regions we want coverage for (whole contigs >= 400bp, excluding 150bp at each end)
cat seq_info.txt | awk '{if ($2 >= 400) {print $1 "\t150\t" $2-150}}' > regions.txt

sambamba depth region -t 8 -L regions.txt ${SPECIES}.bam > depth_info.txt

###############################
#Count total amount of sequence data

cat raw_data/all_F.fastq | awk 'BEGIN {t=0} {if ((NR % 4) == 2) {t += length($1)}} END {print t}' > size_F.txt
cat raw_data/all_R.fastq | awk 'BEGIN {t=0} {if ((NR % 4) == 2) {t += length($1)}} END {print t}' > size_R.txt
cat size_F.txt size_R.txt | awk 'BEGIN {t=0} {print $1; t += $1} END {print t}' > full_size.txt



