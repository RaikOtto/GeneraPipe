GSEA.NormalizeCols <- function(V) { 
  #
  # Stardardize columns of a gene expression matrix
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
  
  col.mean <- apply(V, MARGIN=2, FUN=mean)
  col.sd <- apply(V, MARGIN=2, FUN=sd)
  col.n <- length(V[1,])
  for (i in 1:col.n) {
    if (col.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[,i] <- (V[,i] - col.mean[i])/col.sd[i]
    }
  }
  return(V)
}