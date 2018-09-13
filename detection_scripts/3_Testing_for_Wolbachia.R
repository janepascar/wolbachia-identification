
is.positive <- function(id, pid.threshold=95, length.threshold=98, count.threshold=3) {
  my.data <- read.table(file=paste(id, ".csv", sep=""), col.names =c("read", "ref.id", "pid", "qstart", "qstop", "rstart", "rstop","score", "totalscore"), comment.char="", skip=3)
  my.data$length <- abs(my.data$qstop - my.data$qstart) + 1
  good.alignments <- sum((my.data$length >= length.threshold) & (my.data$pid >= pid.threshold))
  if (good.alignments >= count.threshold) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

args <- commandArgs(trailingOnly = TRUE)
inputfile <- args[1]
outputfile <- args[2]
print(inputfile)
print(outputfile)
acc.list <- read.table(file=inputfile)
acc.list$results <- NA
for (r in 1:nrow(acc.list)) {
  cur.acc <- acc.list[r,1]
  print(cur.acc)
  cur.file <- paste(cur.acc, ".csv", sep="")
  if (file.exists(cur.file)) {
  	acc.list$results[r] <- is.positive(cur.acc)
  }
}

print(acc.list)


write.table(acc.list, file = outputfile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
