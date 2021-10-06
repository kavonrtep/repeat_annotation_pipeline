#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
g = import(commandArgs(T)[1])
gd = mcols(g)
col_name=commandArgs(T)[2]
col_value=commandArgs(T)[3]
inc = gd[,col_name] %in% col_value
g_part = sort(sortSeqlevels(g[inc]))
export(g_part, format = 'gff3', commandArgs(T)[4])
