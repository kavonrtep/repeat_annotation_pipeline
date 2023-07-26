#!/usr/bin/env Rscript
## input 1 - fasta with CLXContig names
## input 2 - annotation
## output 3 - annotated fasta
suppressPackageStartupMessages(library(Biostrings))

clean_contigs = function(s) {
  ## remove all N
  sc = as.character(s)
  sc_trimmed = gsub("N+$", "", gsub("^N+", "", s))
  ## remove blank and short sequences:
  sc_trimmed_not_empty = sc_trimmed[nchar(sc_trimmed) != 0]
  sc_trimmed_short = sc_trimmed_not_empty[nchar(sc_trimmed_not_empty) <= 20]
  sc_trimmed_long = sc_trimmed_not_empty[nchar(sc_trimmed_not_empty) > 20]
  sc_trimmed_short_tarean = sc_trimmed_short[grep("sc_", names(sc_trimmed_short), fixed = TRUE)]
  sc_out = DNAStringSet(c(sc_trimmed_long, sc_trimmed_short_tarean))
}

## annotate_rm_fasta.R input.fasta annot.csv output.fasta
## input 1 - input.fasta - contigs from clustering
## input 2 - CLUSTER_TABLE.csv
##
## output - clean conntigs with appended annotation

## find header row of annotation table
x = readLines(commandArgs(T)[2])

## TODO - check mandatory names!!!
hl = intersect(grep("cluster", tolower(x)), grep("automatic_annotation", tolower(x)))
message("using line ", hl, " as header")

annotation_table = read.table(commandArgs(T)[2], sep = "\t", header = TRUE, skip = hl - 1)
colnames(annotation_table) = tolower(colnames(annotation_table))

contigs = readDNAStringSet(commandArgs(T)[1])


if ("final_annotation" %in% colnames(annotation_table) & all(!is.na(annotation_table$final_annotation))) {
  annot_dict = annotation_table$final_annotation
  message("using final annotation column")
}else {
  message("using automatic annotation column")
  annot_dict = annotation_table$automatic_annotation
}


names(annot_dict) = paste0("CL", annotation_table$cluster)
#print(annot_dict)

contigs_ok = clean_contigs(contigs)
contig_name = gsub("Contig.+", "", names(contigs_ok))

## keep only contigs which are in annot table
include = contig_name %in% names(annot_dict)

contig_name_inc = contig_name[include]
contig_ok_include = contigs_ok[include]

new_name_with_annot = paste0(names(contig_ok_include), "#", annot_dict[contig_name_inc])
names(contig_ok_include) = new_name_with_annot

writeXStringSet(contig_ok_include, filepath = commandArgs(T)[3])
