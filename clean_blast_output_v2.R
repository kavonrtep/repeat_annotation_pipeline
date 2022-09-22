#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(parallel))

convert_names <- function(n, old_sep = "|" , new_sep = "\""){
  # remove all characters which are new_sep with -
  n_new = gsub(old_sep, new_sep,
               gsub(new_sep,"-", n, fixed = TRUE),
               fixed = TRUE)
  return(n_new)
}


infile <-  commandArgs(T)[1]
outfile <- commandArgs(T)[2]

if (FALSE){
    infile = "./tmp/test_output16/blastn_rm_like.tsv"
    outfile = "./tmp/test_out.gff3"
}

# the order of column in blast is:
# fmt2="6 bitscore pident gapopen evalue qseqid qstart qend slen sstrand qcovs sseqid sstart send"

rm_out = read.table(infile, as.is=TRUE, sep="\t", header=FALSE, comment.char="")

gff = GRanges(seqnames = rm_out$V5, ranges = IRanges(start = rm_out$V6, end=rm_out$V7))

# classification sepatator - __ or #
if (grepl("__", rm_out$V11[1], fixed = TRUE)){
  rm_out$V11 <- gsub("__.+","",rm_out$V11)
}else{
  rm_out$V11 <- gsub("^.+#","",rm_out$V11)
}
# repeat class after # symbol - syntax 1
# detect separator
# if "|" is present replace "|" -> "/" and "/" -> "-"
if (any(grepl("|", rm_out$V11, fixed = TRUE))){
  gff$Name <- convert_names(rm_out$V11, old_sep = "|", new_sep = "/")
  message('replacing classification separator character "|" with "/"')
  print(gff)
}else{
  gff$Name <- rm_out$V11
}


gff$bitscore <- rm_out$V1
gff$pident <- rm_out$V2
gff$gapopen <- rm_out$V3
gff$evalue <- rm_out$V4
gff$slen <- rm_out$V8

# merge overlapping regions
merge_neighbors <- function(g){
  rg <- reduce(g, with.revmap=TRUE)
  max_pident <- sapply(rg$revmap, FUN = function(x)max(unlist(g$pident[x])))
  mean_pident <- sapply(rg$revmap, FUN = function(x)weighted.mean(sapply(g$pident[x],max), width(g[x])))
  rg$max_pident <- max_pident
  rg$mean_pident <- mean_pident
  rg$revmap <- NULL
  rg$Name <- g$Name[1]
  rg
}

gff_parts <- split(gff, gff$Name)

result <- unlist(GRangesList(mclapply(gff_parts, merge_neighbors, mc.cores = 8, mc.preschedule = FALSE)))



# split into all overlaping regions
gff_dis <- disjoin(result, with.revmap=TRUE)

gff_dis$score <- mclapply(as.list(gff_dis$revmap), FUN = function(x)gff$bitscore[x], mc.cores = 8, mc.preschedule = FALSE)
gff_dis$pident <- mclapply(as.list(gff_dis$revmap), FUN = function(x)gff$pident[x], mc.cores = 8, mc.preschedule = FALSE)
print(2)
gff_dis$Name <- mclapply(as.list(gff_dis$revmap), FUN = function(x)gff$Name[x], mc.cores = 8)

Final_Classification <- character(length(gff_dis))
# unique names in region
print(3)
N <- lapply(gff_dis$Name, unique)
# number of disctinct classification in each region
L <- sapply(N, length)
Final_Classification[L == 1] <- unlist(N[L==1])
print(4)

gff_dis_unresolved <- gff_dis[L>1]

print(gff_dis_unresolved)
print('lengths')
print(length(gff_dis_unresolved))
print(length(gff_dis_unresolved$Name))
print(length(gff_dis_unresolved$score))

# select classification with best name
best_N <- mapply(x = gff_dis_unresolved$Name, y = gff_dis_unresolved$score, FUN = function(x,y){
  x[which.max(y)]
})
print('-----------------')
head(best_N)
length(Final_Classification)
length(best_N)
print("best_N")
head(best_N)
table(L > 1)


Final_Classification[L>1] <- unlist(best_N)

print('xx-----------------')
gff_dis$Name <- Final_Classification
gff_dis_out <- gff_dis

gff_dis_out$revmap <- NULL
gff_dis_out$score <- sapply(gff_dis_out$score, paste, collapse = ",")
gff_dis_out$pidednt <- sapply(gff_dis_out$pident, paste, collapse = ",")

print(gff_dis_out)
export(gff_dis_out, con=paste0(outfile,"_disjoint.gff"), format = "gff3")

print('xxxxx')


merge_neighbors <- function(g){
  rg <- reduce(g, with.revmap=TRUE)
  max_pident <- sapply(rg$revmap, FUN = function(x)max(unlist(g$pident[x])))
  mean_pident <- sapply(rg$revmap, FUN = function(x)weighted.mean(sapply(g$pident[x],max), width(g[x])))
  rg$max_pident <- max_pident
  rg$mean_pident <- mean_pident
  rg$revmap <- NULL
  rg$Name <- g$Name[1]
  rg
}

gff_dis_parts <- split(gff_dis, Final_Classification)

result <- unlist(GRangesList(mclapply(gff_dis_parts, merge_neighbors, mc.cores = 8, mc.preschedule = FALSE)))

result$ID <- NA
result$revmap <- NA
print(result)
export(result, outfile, format="gff3")



