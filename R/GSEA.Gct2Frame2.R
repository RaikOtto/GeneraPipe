GSEA.Gct2Frame2 <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
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
  content <- readLines(filename)
  content <- content[-1]
  content <- content[-1]
  col.names <- noquote(unlist(strsplit(content[1], "\t")))
  col.names <- col.names[c(-1, -2)]
  num.cols <- length(col.names)
  content <- content[-1]
  num.lines <- length(content)
  
  
  row.nam <- vector(length=num.lines, mode="character")
  row.des <- vector(length=num.lines, mode="character")
  m <- matrix(0, nrow=num.lines, ncol=num.cols)
  
  for (i in 1:num.lines) {
    line.list <- noquote(unlist(strsplit(content[i], "\t")))
    row.nam[i] <- noquote(line.list[1])
    row.des[i] <- noquote(line.list[2])
    line.list <- line.list[c(-1, -2)]
    for (j in 1:length(line.list)) {
      m[i, j] <- as.numeric(line.list[j])
    }
  }
  ds <- data.frame(m)
  names(ds) <- col.names
  row.names(ds) <- row.nam
  return(ds)
}