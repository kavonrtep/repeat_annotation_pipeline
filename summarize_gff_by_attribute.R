#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
g = import(commandArgs(T)[1])
attribute_name = commandArgs(T)[2]

m = mcols(g)
w = width(g)
total_lengths = by(w, INDICES=m[,attribute_name] , sum)
total_counts =  by(w, INDICES=m[,attribute_name] , length)
d = data.frame(attribute = names(total_counts), cbind(counts = total_counts, length=total_lengths))
colnames(d)[1] = attribute_name
d = d[order(d$length, decreasing = TRUE),]
write.table(d, sep = "\t", row.names = FALSE, quote = FALSE)
