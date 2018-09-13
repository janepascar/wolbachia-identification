#PBS -k o
#PBS -l nodes=1:ppn=8
#PBS -l walltime=3:00:00
#PBS -N 6b_contig_selection
#PBS -j oe

# This script will combine the assemblies for each kmer value for each gene and sample
# Once the assemblies are combined in to one file they can be compared to see which kmer is best for each individual sample

# ftsZ
# Make the list of positive accession numbers a variable
SAMPLE=~/downloads/positive_samples_subset.txt

while read SAMPLE; do
  echo ${SAMPLE}
  for kmer in 21 31 41 51
  do
    # This will be the output from the command
    CUR_FILE=~/assembly_contig_selection/combined_${SAMPLE}_ftsZ.fasta
    touch ${CUR_FILE}
    # Make the assemblies a variable 
    CUR_INPUT_FILE=~/assembly_velvet/ftsZ/${SAMPLE}_ftsZ_${kmer}/contigs.fa
    # This will name the contigs so we can know which assembly they came from
    sed "s/>/>${SAMPLE}_k_${kmer}_/" ${CUR_INPUT_FILE} >> ${CUR_FILE}
  done    
done <${SAMPLE}

# groE
SAMPLE=~/downloads/positive_samples_subset.txt

while read SAMPLE; do
  echo ${SAMPLE}
  for kmer in 21 31 41 51
  do
    CUR_FILE=~/assembly_contig_selection/combined_${SAMPLE}_groE.fasta
    touch ${CUR_FILE}
    CUR_INPUT_FILE=~/assembly_velvet/groE/${SAMPLE}_groE_${kmer}/contigs.fa
    sed "s/>/>${SAMPLE}_k_${kmer}_/" ${CUR_INPUT_FILE} >> ${CUR_FILE}
  done    
done <${SAMPLE}

# wsp  
SAMPLE=~/downloads/positive_samples_subset.txt

while read SAMPLE; do
  echo ${SAMPLE}
  for kmer in 21 31 41 51
  do
    CUR_FILE=~/assembly_contig_selection/combined_${SAMPLE}_wsp.fasta
    touch ${CUR_FILE}
    CUR_INPUT_FILE=~/assembly_velvet/wsp/${SAMPLE}_wsp_${kmer}/contigs.fa
    sed "s/>/>${SAMPLE}_k_${kmer}_/" ${CUR_INPUT_FILE} >> ${CUR_FILE}
  done    
done <${SAMPLE}

echo ' done with contig selection ****************************************************************************************************************'

# To find the best contig for each gene they will be BLAST against a database
# The database is specific to the gene being looked at
# For example: ftsZ contigs will be BLAST against an ftsZ database, etc. 
# blast ftsZ 
module load BLAST2.28
SAMPLE=~/downloads/positive_samples_subset.txt
makeblastdb -in ~/sample_genes/ftsZ_genes.fa -dbtype nucl -out ~/sample_genes/ftsZ_DB
while read SAMPLE; do
  echo ${SAMPLE}
  blastn -query ~/assembly_contig_selection/combined_${SAMPLE}_ftsZ.fasta -db ~/sample_genes/ftsZ_DB -outfmt '6 qseqid sseqid length bitscore evalue qstart qend' > ~/contig_BLAST_results/${SAMPLE}_ftsZ.blast
done <${SAMPLE}

# blast groE
module load BLAST2.28
SAMPLE=~/downloads/positive_samples_subset.txt
makeblastdb -in ~/sample_genes/groE_genes.fa -dbtype nucl -out ~/sample_genes/groE_DB
while read SAMPLE; do
  echo ${SAMPLE}
  blastn -query ~/assembly_contig_selection/combined_${SAMPLE}_groE.fasta -db ~/sample_genes/groE_DB -outfmt '6 qseqid sseqid length bitscore evalue qstart qend' > ~/contig_BLAST_results/${SAMPLE}_groE.blast
done <${SAMPLE}

# blast wsp
module load BLAST2.28
SAMPLE=~/downloads/positive_samples_subset.txt
makeblastdb -in ~/sample_genes/wsp_genes.fa -dbtype nucl -out ~/sample_genes/wsp_DB
while read SAMPLE; do
  echo ${SAMPLE}
  blastn -query ~/assembly_contig_selection/combined_${SAMPLE}_wsp.fasta -db ~/sample_genes/wsp_DB -outfmt '6 qseqid sseqid length bitscore evalue qstart qend' > ~/contig_BLAST_results/${SAMPLE}_wsp.blast
done <${SAMPLE}

echo ' done with BLAST ****************************************************************************************************************'