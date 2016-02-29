GSEA.ReadClsFile <- function(file = "NULL") { 
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length=l, mode="character")
  phen.label <- vector(length=l, mode="numeric")
  class.v <- vector(length=s, mode="numeric")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
    phen.label[i] <- i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {
        class.v[i] <- phen.label[j]
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}