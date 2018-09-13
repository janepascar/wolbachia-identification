#!/bin/bash
#PBS -k o 
#PBS -l nodes=1:ppn=8,vmem=32gb,walltime=36:00:00
#PBS -M cholden23@gmail.com
#PBS -m abe 
#PBS -N XXSPECIESXX
#PBS -j oe 


#####################################################
#####################################################

SPECIES=XXSPECIESXX
WORKDIR=/N/dc2/scratch/chrichan/wolbachia/${SPECIES}
ACCESSIONS='XXACCESSIONXX'
MIRAK=31
MINBLASTPIDENT=70
MINBLASTLENGTH=100
MINBLASTEVALUE=1e-10


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
#Download Wolbachia genomes
mkdir -p wo_ref
cd wo_ref

#wPip = group B from Culex quinquefasciatus
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/073/005/GCF_000073005.1_ASM7300v1/GCF_000073005.1_ASM7300v1_genomic.fna.gz
gunzip GCF_000073005.1_ASM7300v1_genomic.fna.gz
mv GCF_000073005.1_ASM7300v1_genomic.fna wPip.fasta

#wMel = group A from Drosophila melanogaster
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/025/GCF_000008025.1_ASM802v1/GCF_000008025.1_ASM802v1_genomic.fna.gz
gunzip GCF_000008025.1_ASM802v1_genomic.fna.gz
mv GCF_000008025.1_ASM802v1_genomic.fna wMel.fasta

#wSim = group ?? from Drosophila simulans
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/376/585/GCF_000376585.1_ASM37658v1/GCF_000376585.1_ASM37658v1_genomic.fna.gz
gunzip GCF_000376585.1_ASM37658v1_genomic.fna.gz
mv GCF_000376585.1_ASM37658v1_genomic.fna wSim.fasta

#wRi = group ?? from Drosophila simulans
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/022/285/GCA_000022285.1_ASM2228v1/GCA_000022285.1_ASM2228v1_genomic.fna.gz
gunzip GCA_000022285.1_ASM2228v1_genomic.fna.gz
mv GCA_000022285.1_ASM2228v1_genomic.fna wRi.fasta

#wOnch = group ?? from Onchocerca volvulus
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/530/755/GCF_000530755.1_W_O_volvulus_Cameroon_v3/GCF_000530755.1_W_O_volvulus_Cameroon_v3_genomic.fna.gz
gunzip GCF_000530755.1_W_O_volvulus_Cameroon_v3_genomic.fna.gz
mv GCF_000530755.1_W_O_volvulus_Cameroon_v3_genomic.fna wOnch.fasta

#wCim = group ?? from Cimex lectularius
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/829/315/GCF_000829315.1_ASM82931v1/GCF_000829315.1_ASM82931v1_genomic.fna.gz
gunzip GCF_000829315.1_ASM82931v1_genomic.fna.gz
mv GCF_000829315.1_ASM82931v1_genomic.fna wCim.fasta

#wTri = group ?? from Trichogramma pretiosum
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/439/985/GCF_001439985.1_wTPRE_1.0/GCF_001439985.1_wTPRE_1.0_genomic.fna.gz
gunzip GCF_001439985.1_wTPRE_1.0_genomic.fna.gz
mv GCF_001439985.1_wTPRE_1.0_genomic.fna wTri.fasta

#Combine them

cat wPip.fasta \
    wMel.fasta \
    wSim.fasta \
    wRi.fasta \
    wOnch.fasta \
    wCim.fasta \
    wTri.fasta > wo_all.fasta

cd ..

###############################
#Pull out Wolbachia-like sequences using mirabait

mkdir -p round1
cd round1

#Initial baiting
mirabait -f fastq -t fastq -k ${MIRAK} ../wo_ref/wo_all.fasta ../raw_data/all_F.fastq initial_bait_F
mirabait -f fastq -t fastq -k ${MIRAK} ../wo_ref/wo_all.fasta ../raw_data/all_R.fastq initial_bait_R

#Make sure to get 'mates' if only one of the members of the pair is a match
cat initial_bait_F.fastq | awk '{if ((NR % 4) == 1) {print substr($1, 2)}}' | sort > F_read_ids.txt
cat initial_bait_R.fastq | awk '{if ((NR % 4) == 1) {print substr($1, 2)}}' | sort > R_read_ids.txt
cat F_read_ids.txt R_read_ids.txt | sort | uniq > all_baited_read_ids.txt
seqtk subseq ../raw_data/all_F.fastq all_baited_read_ids.txt > cur_wo_reads_F.fastq
seqtk subseq ../raw_data/all_R.fastq all_baited_read_ids.txt > cur_wo_reads_R.fastq

cd ..

###############################
#Assemble them

cd round1

module load spades/intel/3.11.1
spades.py -1 cur_wo_reads_F.fastq -2 cur_wo_reads_R.fastq --only-assembler -t 8 -m 30 -o ./wo_assembly

cd ..

###############################
#BLAST against reference Wolbachia genomes and keep only the assembled contigs/scaffolds that have decent hits

cd round1

#Do the blast search
makeblastdb -in ../wo_ref/wo_all.fasta -dbtype nucl
blastn -task dc-megablast -num_threads 8 -query wo_assembly/scaffolds.fasta -db ../wo_ref/wo_all.fasta -outfmt '6 qseqid sseqid evalue bitscore pident length' -evalue 1e-10 > scaffold_blast_results.txt

#Filter down to reliable BLAST hits
cat scaffold_blast_results.txt | awk -v minevalue=${MINBLASTEVALUE} -v minlen=${MINBLASTLENGTH} -v minpident=${MINBLASTPIDENT} '{if (($3 <= minevalue) && ($5 >= minpident) && ($6 >= minlen)) {print $1}}' | sort | uniq > wo_matching_scaffolds.txt

#Get a fasta file with only the matching contigs
seqtk subseq wo_assembly/scaffolds.fasta wo_matching_scaffolds.txt > wo_matching_scaffolds.fasta

cd ..

###############################
#Repeat the mapping, read extraction, & assembly process a few more times

PREVROUND=1
for CURROUND in 2 3 4 5
do
    mkdir -p round${CURROUND}
    cd round${CURROUND}

    #Extract reads
    PREVASSEMBLY=../round${PREVROUND}/wo_matching_scaffolds.fasta
    mirabait -f fastq -t fastq -k ${MIRAK} ${PREVASSEMBLY} ../raw_data/all_F.fastq initial_bait_F
    mirabait -f fastq -t fastq -k ${MIRAK} ${PREVASSEMBLY} ../raw_data/all_R.fastq initial_bait_R
    cat initial_bait_F.fastq | awk '{if ((NR % 4) == 1) {print substr($1, 2)}}' | sort > F_read_ids.txt
    cat initial_bait_R.fastq | awk '{if ((NR % 4) == 1) {print substr($1, 2)}}' | sort > R_read_ids.txt
    cat F_read_ids.txt R_read_ids.txt | sort | uniq > all_baited_read_ids.txt
    seqtk subseq ../raw_data/all_F.fastq all_baited_read_ids.txt > cur_wo_reads_F.fastq
    seqtk subseq ../raw_data/all_R.fastq all_baited_read_ids.txt > cur_wo_reads_R.fastq
    
    #Assemble them
    spades.py -1 cur_wo_reads_F.fastq -2 cur_wo_reads_R.fastq --only-assembler -t 8 -m 30 -o ./wo_assembly
    
    #Blast and keep only the ones that show some similarity to other Wolbachia genomes
    blastn -task dc-megablast -num_threads 8 -query wo_assembly/scaffolds.fasta -db ../wo_ref/wo_all.fasta -outfmt '6 qseqid sseqid evalue bitscore pident length' -evalue 1e-10 > scaffold_blast_results.txt
    cat scaffold_blast_results.txt | awk -v minevalue=${MINBLASTEVALUE} -v minlen=${MINBLASTLENGTH} -v minpident=${MINBLASTPIDENT} '{if (($3 <= minevalue) && ($5 >= minpident) && ($6 >= minlen)) {print $1}}' | sort | uniq > wo_matching_scaffolds.txt
    seqtk subseq wo_assembly/scaffolds.fasta wo_matching_scaffolds.txt > wo_matching_scaffolds.fasta
    
    #Get set up for next round
    PREVROUND=${CURROUND}
    
    cd ..
done

###############################
#Evaluate the assemblies 

module load python/2.7.13
QUAST=/N/dc2/scratch/chrichan/temp/quast-4.5/quast.py

for CUR_REF_BASE in wMel wPip
do
    CUR_WO_REF=wo_ref/${CUR_REF_BASE}.fasta
    ${QUAST} -R ${CUR_WO_REF} -t 8 \
      -o ${CUR_REF_BASE}_eval \
      round1/wo_matching_scaffolds.fasta \
      round2/wo_matching_scaffolds.fasta \
      round3/wo_matching_scaffolds.fasta \
      round4/wo_matching_scaffolds.fasta \
      round5/wo_matching_scaffolds.fasta
done


###############################
#Generate histograms

echo 'all.data <- read.table("cov_data.txt", col.names=c("id", "cov", "length"), header=F)' > plot_cov.R
echo 'tot.len <- sum(all.data$length)' >> plot_cov.R
echo 'pdf(file="coverage_histogram.pdf", width=5, height=3)' >> plot_cov.R
echo 'par(mar=c(4,4,1,1))' >> plot_cov.R
echo 'n.bins <- 300' >> plot_cov.R
echo 'n.axis.skips <- 3' >> plot_cov.R
echo 'proportion.plot <- 0.98' >> plot_cov.R
echo 'mybins <- cut(all.data$cov, n.bins)' >> plot_cov.R
echo 'plot.data <- tapply(X=all.data$length, INDEX=mybins, FUN=sum)' >> plot_cov.R
echo 'plot.data <- ifelse(is.na(plot.data), 0, plot.data)' >> plot_cov.R
echo 'cumsum.data <- cumsum(plot.data)' >> plot_cov.R
echo 'xaxis.max <- round(sum(plot.data) * proportion.plot)' >> plot_cov.R
echo 'plot.data <- plot.data[cumsum.data <= xaxis.max]' >> plot_cov.R
echo 'barplot(height=plot.data/1000, space=0, border=NA, xlab="Coverage", ylab="Span (kb)", xaxt="n", main=sprintf("Total span = %.3f Mb", round(tot.len/1e6, digits=3)))' >> plot_cov.R
echo 'bin.labels <- gsub(pattern="\\]", replacement="", names(plot.data))' >> plot_cov.R
echo 'bin.labels <- gsub(pattern="\\(", replacement="", bin.labels)' >> plot_cov.R
echo 'bin.means <- round(unlist(lapply(strsplit(bin.labels, split=","), FUN=function(X){mean(as.numeric(X))})), digits=1)' >> plot_cov.R
echo 'mybins.show <- seq(from=1, to=length(bin.means), by=n.axis.skips)' >> plot_cov.R
echo 'axis(side=1, at=(1:nrow(plot.data))[mybins.show], labels=bin.means[mybins.show], cex.axis=0.8, las=2)' >> plot_cov.R
echo 'dev.off()' >> plot_cov.R

module load r/3.3.1

for CURROUND in 1 2 3 4 5
do
    cd round${CURROUND}
    grep '^>' wo_matching_scaffolds.fasta | sed 's/>//' | awk -F '_' '{print $0 "\t" $6 "\t" $4}' > cov_data.txt
    Rscript --vanilla ../plot_cov.R
    cd ..
done


###############################
#Count total amount of sequence data

cat raw_data/all_F.fastq | awk 'BEGIN {t=0} {if ((NR % 4) == 2) {t += length($1)}} END {print t}' > size_F.txt
cat raw_data/all_R.fastq | awk 'BEGIN {t=0} {if ((NR % 4) == 2) {t += length($1)}} END {print t}' > size_R.txt
cat size_F.txt size_R.txt | awk 'BEGIN {t=0} {print $1; t += $1} END {print t}' > full_size.txt



