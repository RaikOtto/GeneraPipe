GSEA.NormalizeRows <- function(V) { 
  #
  # Stardardize rows of a gene expression matrix
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
  
  row.mean <- apply(V, MARGIN=1, FUN=mean)
  row.sd <- apply(V, MARGIN=1, FUN=sd)
  row.n <- length(V[,1])
  for (i in 1:row.n) {
    if (row.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[i,] <- (V[i,] - row.mean[i])/row.sd[i]
    }
  }
  return(V)
}