#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
gff = import(commandArgs(T)[1])
BN = gsub("[.]gff(3){0,1}","",basename(commandArgs(T)[1]))
min_width = as.numeric(commandArgs(T)[2])
outdir = commandArgs(T)[3]



gff_min_width=gff[width(gff)>=min_width]



dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

gff_min_width_parts = split(gff_min_width, f=gff_min_width$Name)

x = sapply(names(gff_min_width_parts),
       function(x) export(
                     gff_min_width_parts[[x]],
                     format="gff3",
                     con=paste0(
                       outdir,"/",
                       gsub("/","_",x),
                       "_", min_width ,"plus.gff3"
                     )
                   ))


