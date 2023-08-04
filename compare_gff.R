#!/usr/bin/env Rscript
# function from questionr R package:
wtable <- function(x, y = NULL, weights = NULL, digits = 3, normwt = FALSE, useNA = c("no", "ifany", "always"), na.rm = TRUE, na.show = FALSE, exclude = NULL)
{
  if (is.null(weights)) {
    warning("no weights argument given, using uniform weights of 1")
    weights <- rep(1, length(x))
  }
  if (length(x) != length(weights)) stop("x and weights lengths must be the same")
  if (!is.null(y) & (length(x) != length(y))) stop("x and y lengths must be the same")
  miss.usena <- missing(useNA)
  useNA <- match.arg(useNA)
  weights[is.na(weights)] <- 0
  if (normwt) {
    weights <- weights * length(x)/sum(weights)
  }

  if (!missing(na.show) || !missing(na.rm)) {
    warning("'na.rm' and 'na.show' are ignored when 'useNA' is provided.")
  }
  if (useNA != "no" || (na.show && miss.usena)) {
    if (match(NA, exclude, nomatch = 0L)) {
      warning("'exclude' containing NA and 'useNA' != \"no\"' are a bit contradicting")
    }
    x <- addNA(x)
    if (!is.null(y)) y <- addNA(y)
  }
  if (useNA == "no" || (na.rm && miss.usena)) {
     s <- !is.na(x) & !is.na(weights)
     if (!is.null(y)) s <- s & !is.na(y)
     x <- x[s, drop = FALSE]
     if (!is.null(y)) y <- y[s, drop = FALSE]
     weights <- weights[s]
  }
  if (!is.null(exclude)) {
    s <- !(x %in% exclude)
    if (!is.null(y)) s <- s & !(y %in% exclude)
    x <- factor(x[s, drop = FALSE])
    if (!is.null(y)) y <- factor(y[s, drop = FALSE])
    weights <- weights[s]
  }
  if (is.null(y)) {
    result <- tapply(weights, x, sum, simplify = TRUE)
  }
  else {
    result <- tapply(weights, list(x,y), sum, simplify = TRUE)
  }
  result[is.na(result)] <- 0
  tab <- as.table(result)
  if (useNA == "ifany") {
    if (!is.null(y)) {
      if (sum(tab[,is.na(colnames(tab))]) == 0) tab <- tab[,!is.na(colnames(tab))]
      if (sum(tab[is.na(rownames(tab)),]) == 0) tab <- tab[!is.na(rownames(tab)),]
    } else {
      if (tab[is.na(names(tab))] == 0) tab <- tab[!is.na(names(tab))]
    }
  }
  tab
}

suppressPackageStartupMessages(library(rtracklayer))
library(parallel)
g1 <- import(commandArgs(T)[1], format = "GFF")
g2 <-  import(commandArgs(T)[2], format = "GFF")
# assume non inrarange overlaps!
attribute_name <- commandArgs(T)[3]
if (FALSE){
  g1 <- import("test_data/RM_LTR_TE.gff3")
  g2 <- import("test_data/RM_RE.gff3")
  g2 <- import("test_data/RM_RE_CL.gff3")
  g2 <- import("test_data/RM_RE_v3.gff3")
  g2 <- import("test_data/RM_RE_v4.gff3")
  g1 <- import("test_data/RM_RE_v4.gff3")

  attribute_name <- "Name"
}

g1$SOURCE <- 1
g2$SOURCE <- 2

g12 <-  append(g1,g2)
g12_disjoin <- disjoin(g12, ignore.strand=TRUE, with.revmap=TRUE)

g12$sel_attr <- mcols(g12)[,attribute_name]

annot_name <- mclapply(as.list(g12_disjoin$revmap), FUN = function(x)g12$sel_attr[x], mc.cores = 8)
annot_source <- mclapply(as.list(g12_disjoin$revmap), FUN = function(x)g12$SOURCE[x], mc.cores = 8)

c1 <- sapply(annot_source, function (x) 1 %in% x & !(2 %in% x))
c2 <- sapply(annot_source, function (x) 2 %in% x & !(1 %in% x))
c12 <- !(c1 | c2)


n1 <- cbind(sapply(annot_name[c1], "[", 1 ), "No Annotation")
n2 <- cbind("No Annotation", sapply(annot_name[c2], "[", 1))

n12 <- do.call(rbind, annot_name[c12])

n12all <- matrix(character(), ncol = 2, nrow = length(g12_disjoin))

if (any(c1)){
  n12all[c1,] <- n1
}
if (any(c2)){
  n12all[c2,] <- n2
}
if (any(c12)){
  n12all[c12,] <- n12
}

traw <- wtable(n12all[,1], n12all[,2], weights = width(g12_disjoin))
tbl_df <- as.data.frame((traw))
tbl_df <- tbl_df[order(tbl_df$Freq, decreasing = TRUE),]
colnames(tbl_df) <-  c("Annotation1", "Annotation2", "Overlap[bp]")

## NOTE - if name change - modify it in xml too
write.table(tbl_df, file = 'annotation_overlap_long.csv', sep="\t", quote=FALSE, row.names = FALSE)

tbl <- as.matrix(traw)
tbl <- tbl[order(rowSums(tbl), decreasing = TRUE),]
tbl <- tbl[, order(colSums(tbl), decreasing = TRUE)]
write.table(tbl, file= 'annotation_overlap.csv', sep="\t", quote = FALSE)

