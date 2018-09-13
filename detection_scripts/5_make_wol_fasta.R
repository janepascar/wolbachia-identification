args <- commandArgs(trailingOnly = TRUE)

dict.file <- args[1]
input.file <- args[2]
output.base <- args[3]

sam.colnames <- c("query", "flag", "ref", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual", "nh", "as", "nm")

sam.data <- read.table(input.file, col.names=sam.colnames, comment.char="@", stringsAsFactors=F)

#Function to parse CIGAR strings to get the total number of M's (number of aligned bases)
get.num.m <- function(cigar.string) {
	nums <- strsplit(x=cigar.string, split="[SMHID]", fixed=F)[[1]]
	letters <- strsplit(x=cigar.string, split="[0-9]+", fixed=F)[[1]]
	letters <- letters[nchar(letters) >= 1]
	m.vector <- as.numeric(nums[letters=="M"])
	return(sum(m.vector))
}

sam.data$num.m <- sapply(X=sam.data$cigar, FUN=get.num.m)

#Sort the data
sam.data <- sam.data[order(sam.data$num.m, decreasing=T),]
sam.data <- sam.data[order(sam.data$query),]


#PROBLEM: In the sequence column of the SAM file, MagicBlast just outputs the forward and reverse reads concatenated
#         instead of treating them as separate reads!
#If we just take this sequence and use it in our assembly, things will get screwed up
#So, let's figure out if this dataset looks like it consisted of single- or double-ended reads
#If it was double-ended, we'll split the sequences in the middle when we output them

#How long is the longest read in this dataset?
read.len <- max(nchar(sam.data$seq), na.rm=T)

#Only keep the best hit for each read
keep.data <- sam.data[!duplicated(sam.data$query),]

#Only keep reads that have at least a minimum number of matches to the reference
#Figure out how many M's we need to keep a read -- let's not accept anything less than 50, and if we have long reads, let's
# only keep those that are at least 90% of half of the read length
keep.thresh <- max(c(50, floor((read.len/2)*0.9)))
keep.data <- keep.data[keep.data$num.m >= keep.thresh,]

#Identify which gene each read aligns to
my.dict <- read.table(dict.file, col.names=c("num", "gene"), stringsAsFactors=F)
keep.data$gene <- my.dict$gene[match(x=keep.data$ref, table=my.dict$num)]


num.whole.read.matches <- sum(keep.data$num.m >= (read.len * 0.9))

write.fasta <- function(dataset, filename, split.reads=F) {
	if (split.reads) {
		output.vector <- rep("", nrow(dataset)*4)
		split.length <- floor(read.len / 2)
	} else {
		output.vector <- rep("", nrow(dataset)*2)
	}
	n.seqs <- length(output.vector) / 2
	fasta.names <- seq(from=1, to=length(output.vector), by=2)
	fasta.seqs <- fasta.names+1
	output.vector[fasta.names] <- paste(">Seq_", as.character(seq(from=1, to=n.seqs, by=1)), sep="")
	cur.seq <- 1
	for (n in seq(from=1, to=nrow(dataset), by=1)) {
		if (split.reads) {
			cur.read <- dataset$seq[n]
			output.vector[fasta.seqs[cur.seq]] <- substr(cur.read, 1, split.length)
			cur.seq <- cur.seq + 1
			output.vector[fasta.seqs[cur.seq]] <- substr(cur.read, split.length + 1, nchar(cur.read))
			cur.seq <- cur.seq + 1
		} else {
			cur.read <- dataset$seq[n]
			output.vector[fasta.seqs[cur.seq]] <- cur.read
		}
	}
	write(x=output.vector, file=filename)
}

all.genes <- c("wsp", "groE", "ftsZ")

for (cur.gene in all.genes) {
	cur.output.file <- paste(output.base, "_", cur.gene, ".fasta", sep="")
	cur.data <- keep.data[keep.data$gene == cur.gene,]
	
	#If there are at least 20 sequence reads that match along nearly their entire length, then
	# let's assume this dataset is single-ended, and output whole reads
	#If not, we'll split the reads in the middle before we output them
	write.fasta(cur.data, cur.output.file, (num.whole.read.matches < 20))
}

