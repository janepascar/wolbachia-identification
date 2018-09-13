# This script is used to download Sequence Read Archive entries based on accession number to a cluster
# Use the 'screen' command to download data in the background

# Change to the directory where the downloaded data will be stored  
cd ~/downloads

# List the accession numbers in a text file
ACC=~/downloads/accession_numbers.txt

# Start a loop to download all SRA accession numbers listed in the text file
# A maximum of 50000000 spots will be downloaded
# The accession number that is currently being downloaded will be echoed to the screen
# The downloaded sequences will be stored in a new directory in the ~/downloads directory
mkdir SRA_genome_sequences
while read CUR_ACC; do 
  echo ${CUR_ACC}
  ~/bin/fastq-dump.2.8.0 --readids --outdir SRA_genome_sequences --origfmt --skip-technical -X 50000000 ${CUR_ACC}
done <${ACC}

echo 'done ******************************************************************************'

