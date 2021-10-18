#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
g = import(commandArgs(T)[1], format = "GFF")
attribute_name = commandArgs(T)[2]

m = mcols(g)
AN = commandArgs(T)[2]
## check if attribute is in gff:
if (!AN %in% colnames(m)){
  stop(paste0("attribute ", AN," not in input gff\n", "Available attributes are:\n", paste(colnames(m), collapse = "\n ")))
}

w = width(g)
total_lengths = by(w, INDICES=m[,attribute_name] , sum)
total_counts =  by(w, INDICES=m[,attribute_name] , length)
total_length_summary = do.call(rbind,as.list(by(w, INDICES=m[,attribute_name] , summary)))
colnames(total_length_summary)  = paste0("length.", colnames(total_length_summary))

d = data.frame(attribute = names(total_counts), cbind(count = total_counts, total_length=total_lengths, length=total_length_summary))
colnames(d)[1] = attribute_name
d = d[order(d$total_length, decreasing = TRUE),]
write.table(d, sep = "\t", row.names = FALSE, quote = FALSE, file = commandArgs(T)[3])
