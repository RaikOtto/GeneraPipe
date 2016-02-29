GSEA.write.gct <- function (gct, filename) 
{
  f <- file(filename, "w")
  cat("#1.2", "\n", file = f, append = TRUE, sep = "")
  cat(dim(gct)[1], "\t", dim(gct)[2], "\n", file = f, append = TRUE, sep = "")
  cat("Name", "\t", file = f, append = TRUE, sep = "")
  cat("Description", file = f, append = TRUE, sep = "")
  names <- names(gct)
  cat("\t", names[1], file = f, append = TRUE, sep = "")
  for (j in 2:length(names)) {
    cat("\t", names[j], file = f, append = TRUE, sep = "")
  }
  cat("\n", file = f, append = TRUE, sep = "\t")
  oldWarn <- options(warn = -1)
  m <- matrix(nrow = dim(gct)[1], ncol = dim(gct)[2] +  2)
  m[, 1] <- row.names(gct)
  m[, 2] <- row.names(gct)
  index <- 3
  for (i in 1:dim(gct)[2]) {
    m[, index] <- gct[, i]
    index <- index + 1
  }
  write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
  close(f)
  options(warn = 0)
  return(gct)
}