all.data <- read.table("depth_info.txt", col.names=c("id", "start", "end", "reads", "cov", "sample"), header=F)
all.data$length <- all.data$end - all.data$start + 1
tot.len <- sum(all.data$length)
pdf(file="coverage_histogram.pdf", width=5, height=3)
par(mar=c(4,4,1.5,1))
n.bins <- 30
n.axis.skips <- 3
proportion.plot <- 1
all.data$cov.log <- log(all.data$cov)

keep.data <- all.data

#Get weighted median coverage
keep.data.sorted <- keep.data[order(keep.data$cov),]
keep.data.sorted$cumlen <- cumsum(keep.data.sorted$length)
middle.base <- max(keep.data.sorted$cumlen)/2
middle.row <- which(keep.data.sorted$cumlen >= middle.base)[1]
median.cov <- keep.data.sorted$cov[middle.row]

mybins <- cut(keep.data$cov.log, n.bins)
plot.data <- tapply(X=keep.data$length, INDEX=mybins, FUN=sum)
plot.data <- ifelse(is.na(plot.data), 0, plot.data)
cumsum.data <- cumsum(plot.data)
xaxis.max <- round(sum(plot.data) * proportion.plot)
plot.data <- plot.data[cumsum.data <= xaxis.max]
barplot(height=plot.data/1000, space=0, border=NA, xlab="Coverage", ylab="Span (kb)", xaxt="n", main=sprintf("Total span assessed for coverage = %.3f Mb\nMedian coverage = %0.2f", round(tot.len/1e6, digits=3), round(median.cov, digits=2)), cex.main=0.5)
bin.labels <- gsub(pattern="\\]", replacement="", names(plot.data))
bin.labels <- gsub(pattern="\\(", replacement="", bin.labels)
bin.means <- round(unlist(lapply(strsplit(bin.labels, split=","), FUN=function(X){mean(as.numeric(X))})), digits=1)
mybins.show <- seq(from=1, to=length(bin.means), by=n.axis.skips)
bin.means <- exp(bin.means)
pretty.bin.means <- sprintf("%0.1f", round(bin.means, digits=1))
axis(side=1, at=(1:nrow(plot.data))[mybins.show], labels=pretty.bin.means[mybins.show], cex.axis=0.8, las=2)
dev.off()
