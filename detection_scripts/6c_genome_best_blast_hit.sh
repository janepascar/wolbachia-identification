#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -N test_best_blast_hit
#PBS -j oe

module load bedtools

SAMPLE=~/downloads/positive_samples_subset.txt

while read SAMPLE; do
  echo ${SAMPLE}
  for GENE in wsp ftsZ groE
  do
    ORIGINAL_FASTA=~/assembly_contig_selection/combined_${SAMPLE}_${GENE}.fasta
    BLAST_RESULTS=~/contig_BLAST_results/${SAMPLE}_${GENE}.blast

    # Find out what the best BLAST hit is for the current sample and gene
    # NOTE: BLAST sequence coordinates are 1-based; BED files use 0-based
    cat ${BLAST_RESULTS} | awk 'BEGIN {contig=""; best_len=0; start=1; end=1;} {if (best_len < $3) {contig=$1; best_len=$3; start=$6-1; end=$7} } END {print contig "\t" start "\t" end}' > ~/BEST_contigs/${SAMPLE}_${GENE}.bed
    echo " best BLAST hit found for ${SAMPLE}_${GENE} ******************************************************************************"

    # Now extract the aligned part of the sequence from the best BLAST hit
    NEW_FASTA=${SAMPLE}_${GENE}_best.fasta
    bedtools getfasta -fi ${ORIGINAL_FASTA} -bed ~/BEST_contigs/${SAMPLE}_${GENE}.bed -fo ~/BEST_contigs/aligned_sequences/${NEW_FASTA}
  done
done <${SAMPLE}

echo ' done with contig selection ****************************************************************************'

# Put all of the blast hits into one fasta file 
cat *_ftsZ*>>all_ftsZ.fasta
cat *_groE*>>all_groE.fasta
cat *_wsp*>>all_wsp.fasta

mv all* ~/phylogeny_construction/

# Add known positive control samples that will be included in the phylogeny
