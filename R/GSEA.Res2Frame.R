GSEA.Res2Frame <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in RES format and converts it into an R data frame
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
  
  header.cont <- readLines(filename, n = 1)
  temp <- unlist(strsplit(header.cont, "\t"))
  colst <- length(temp)
  header.labels <- temp[seq(3, colst, 2)]
  ds <- read.delim(filename, header=F, row.names = 2, sep="\t", skip=3, blank.lines.skip=T, comment.char="", as.is=T)
  colst <- length(ds[1,])
  cols <- (colst - 1)/2
  rows <- length(ds[,1])
  A <- matrix(nrow=rows - 1, ncol=cols)
  A <- ds[1:rows, seq(2, colst, 2)]
  table1 <- data.frame(A)
  names(table1) <- header.labels
  return(table1)
}