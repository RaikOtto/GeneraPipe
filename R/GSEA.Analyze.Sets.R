GSEA.Analyze.Sets <- function(
  directory,
  topgs = "",
  non.interactive.run = FALSE,
  height = 12,
  width = 17) {
  
  file.list <- list.files(directory)
  files <- file.list[regexpr(pattern = ".report.", file.list) > 1]
  max.sets <- length(files)
  
  set.table <- matrix(nrow = max.sets, ncol = 5)
  
  for (i in 1:max.sets) {
    temp1 <-  strsplit(files[i], split=".report.")
    temp2 <-  strsplit(temp1[[1]][1], split="\\.")
    s <- length(temp2[[1]])
    prefix.name <- paste(temp2[[1]][1:(s-1)], sep="", collapse="")
    set.name <- temp2[[1]][s]
    temp3 <-  strsplit(temp1[[1]][2], split="\\.")
    phenotype <- temp3[[1]][1]
    seq.number <- temp3[[1]][2]
    dataset <- paste(temp2[[1]][1:(s-1)], sep="", collapse=".")
    
    set.table[i, 1] <- files[i]
    
    set.table[i, 3] <- phenotype
    set.table[i, 4] <- as.numeric(seq.number)
    set.table[i, 5] <- dataset
    
    #      set.table[i, 2] <- paste(set.name, dataset, sep ="", collapse="")
    set.table[i, 2] <- paste(i, substr(set.name, 1, 20), sep = "_")
  }
  
  print(c("set name=", prefix.name))
  doc.string <- prefix.name
  
  set.table <- noquote(set.table)
  phen.order <- order(set.table[, 3], decreasing = T)
  set.table <- set.table[phen.order,]
  phen1 <- names(table(set.table[,3]))[1]
  phen2 <- names(table(set.table[,3]))[2]
  set.table.phen1 <- set.table[set.table[,3] == phen1,]
  set.table.phen2 <- set.table[set.table[,3] == phen2,]
  
  seq.order <- order(as.numeric(set.table.phen1[, 4]), decreasing = F)
  set.table.phen1 <- set.table.phen1[seq.order,]
  seq.order <- order(as.numeric(set.table.phen2[, 4]), decreasing = F)
  set.table.phen2 <- set.table.phen2[seq.order,]
  
  #   max.sets.phen1 <- length(set.table.phen1[,1])
  #   max.sets.phen2 <- length(set.table.phen2[,1])
  
  if (topgs == "") {
    max.sets.phen1 <- length(set.table.phen1[,1])
    max.sets.phen2 <- length(set.table.phen2[,1])
  } else {
    max.sets.phen1 <- ifelse(topgs > length(set.table.phen1[,1]), length(set.table.phen1[,1]), topgs) 
    max.sets.phen2 <- ifelse(topgs > length(set.table.phen2[,1]), length(set.table.phen2[,1]), topgs)
  }
  
  # Analysis for phen1
  
  leading.lists <- NULL
  for (i in 1:max.sets.phen1) {
    inputfile <- paste(directory, paste( set.table.phen1[i, 1], sep=""), sep = "/", collapse="")
    gene.set <- read.table(file=inputfile, sep="\t", header=T, comment.char="", as.is=T)
    leading.set <- as.vector(gene.set[gene.set[,"CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set 
    } else{
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow=max.sets.phen1, ncol=max.genes)
  for (i in 1:max.sets.phen1) {
    M[i,] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch = 0))   # notice that the sign is 0 (no tag) or 1 (tag) 
  }
  
  Inter <- matrix(0, nrow=max.sets.phen1, ncol=max.sets.phen1)
  for (i in 1:max.sets.phen1) {
    for (j in 1:max.sets.phen1) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }
  
  Itable <- data.frame(Inter)
  names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen1, sep="" ), sep ="/", collapse="")
      windows(height = width, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  }
  
  GSEA.ConsPlot(Itable, col.names = set.table.phen1[1:max.sets.phen1, 2], main = " ", sub=paste("Leading Subsets Overlap ", doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ")
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  # Save leading subsets in a GCT file
  
  D.phen1 <- data.frame(M)
  names(D.phen1) <- all.leading.genes
  row.names(D.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  output <- paste(directory, paste( doc.string, ".leading.genes.", phen1, ".gct", sep=""), sep = "/" )
  GeneraPipe:::GSEA.write.gct(D.phen1, filename=output)
  
  # Save leading subsets as a single gene set in a .gmt file
  
  row.header <- paste(doc.string, ".all.leading.genes.", phen1, sep="")
  output.line <- paste(all.leading.genes, sep="\t", collapse="\t")
  output.line <- paste(row.header, row.header, output.line, sep="\t", collapse="")
  output <- paste(directory, paste( doc.string, ".all.leading.genes.", phen1, ".gmt", sep=""), sep = "/" )
  write(noquote(output.line), file = output, ncolumns = length(output.line))
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen1, sep=""), sep = "/", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  GeneraPipe:::GSEA.HeatMapPlot2(V = data.matrix(D.phen1), row.names = row.names(D.phen1), col.names = names(D.phen1), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  DT1.phen1 <- data.matrix(t(D.phen1))
  DT2.phen1 <- data.frame(DT1.phen1)
  names(DT2.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
  row.names(DT2.phen1) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen1, filename=outputfile2.phen1)
  
  # Analysis for phen2
  
  leading.lists <- NULL
  for (i in 1:max.sets.phen2) {
    inputfile <- paste(directory, set.table.phen2[i, 1], sep="/", collapse="")
    gene.set <- read.table(file=inputfile, sep="\t", header=T, comment.char="", as.is=T)
    leading.set <- as.vector(gene.set[gene.set[,"CORE_ENRICHMENT"] == "YES", "SYMBOL"])
    leading.lists <- c(leading.lists, list(leading.set))
    if (i == 1) {
      all.leading.genes <- leading.set 
    } else{
      all.leading.genes <- union(all.leading.genes, leading.set)
    }
  }
  max.genes <- length(all.leading.genes)
  M <- matrix(0, nrow=max.sets.phen2, ncol=max.genes)
  for (i in 1:max.sets.phen2) {
    M[i,] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch=0))   # notice that the sign is 0 (no tag) or 1 (tag) 
  }
  
  Inter <- matrix(0, nrow=max.sets.phen2, ncol=max.sets.phen2)
  for (i in 1:max.sets.phen2) {
    for (j in 1:max.sets.phen2) {
      Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], leading.lists[[j]]))
    }
  }
  
  Itable <- data.frame(Inter)
  names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen2, sep=""), sep = "/", collapse="")
      windows(height = width, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.overlap.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = width, width = width)
    }
  }
  
  GeneraPipe:::GSEA.ConsPlot(Itable, col.names = set.table.phen2[1:max.sets.phen2, 2], main = " ", sub=paste("Leading Subsets Overlap ", doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ")
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  # Save leading subsets in a GCT file
  
  D.phen2 <- data.frame(M)
  names(D.phen2) <- all.leading.genes
  row.names(D.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  output <- paste(directory, paste( doc.string, ".leading.genes.", phen2, ".gct", sep=""), sep = "/")
  GeneraPipe:::GSEA.write.gct(D.phen2, filename=output)
  
  # Save primary subsets as a single gene set in a .gmt file
  
  row.header <- paste(doc.string, ".all.leading.genes.", phen2, sep="")
  output.line <- paste(all.leading.genes, sep="\t", collapse="\t")
  output.line <- paste(row.header, row.header, output.line, sep="\t", collapse="")
  output <- paste(directory, paste( doc.string, ".all.leading.genes.", phen2, ".gmt", sep=""), sep = "/" )
  write(noquote(output.line), file = output, ncolumns = length(output.line))
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen2, sep=""),sep ="/", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  GeneraPipe:::GSEA.HeatMapPlot2(V = data.matrix(D.phen2), row.names = row.names(D.phen2), col.names = names(D.phen2), main = "Leading Subsets Assignment", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  DT1.phen2 <- data.matrix(t(D.phen2))
  DT2.phen2 <- data.frame(DT1.phen2)
  names(DT2.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
  row.names(DT2.phen2) <- all.leading.genes
  #   GSEA.write.gct(DT2.phen2, filename=outputfile2.phen2)
  
  # Resort columns and rows for phen1
  
  A <- data.matrix(D.phen1)
  A.row.names <- row.names(D.phen1)
  A.names <- names(D.phen1)
  
  # Max.genes
  
  #   init <- 1
  #   for (k in 1:max.sets.phen1) { 
  #      end <- which.max(cumsum(A[k,]))
  #      if (end - init > 1) {
  #         B <- A[,init:end]
  #         B.names <- A.names[init:end]
  #        dist.matrix <- dist(t(B))
  #         HC <- hclust(dist.matrix, method="average")
  ##         B <- B[,HC$order] + 0.2*(k %% 2)
  #        B <- B[,HC$order] 
  #         A[,init:end] <- B
  #         A.names[init:end] <- B.names[HC$order]
  #         init <- end + 1
  #     }
  #   }
  
  #   windows(width=14, height=10)
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ") 
  
  dist.matrix <- dist(t(A))
  HC <- hclust(dist.matrix, method="average")
  A <- A[, HC$order]
  A.names <- A.names[HC$order]
  
  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method="average")
  A <- A[HC$order,]
  A.row.names <- A.row.names[HC$order]
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen1, sep=""), sep = "/", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  
  
  cmap <-  c("#AAAAFF", "#111166")
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  GeneraPipe:::GSEA.HeatMapPlot2(V = t(A), row.names = A.names, col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen1, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  text.filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen1, ".txt", sep=""), sep = "/", collapse="")
  line.list <- c("Gene", A.row.names)
  line.header <- paste(line.list, collapse="\t")
  line.length <- length(A.row.names) + 1
  write(line.header, file = text.filename, ncolumns = line.length)
  write.table(t(A), file=text.filename, append = TRUE, quote = FALSE, col.names= FALSE, row.names = TRUE, sep = "\t")
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
  # resort columns and rows for phen2
  
  A <- data.matrix(D.phen2)
  A.row.names <- row.names(D.phen2)
  A.names <- names(D.phen2)
  
  # Max.genes
  
  #   init <- 1
  #   for (k in 1:max.sets.phen2) { 
  #      end <- which.max(cumsum(A[k,]))
  #      if (end - init > 1) {
  #         B <- A[,init:end]
  #         B.names <- A.names[init:end]
  #         dist.matrix <- dist(t(B))
  #         HC <- hclust(dist.matrix, method="average")
  ##         B <- B[,HC$order] + 0.2*(k %% 2)
  #        B <- B[,HC$order] 
  #         A[,init:end] <- B
  #         A.names[init:end] <- B.names[HC$order]
  #         init <- end + 1
  #     }
  #   }
  
  #  windows(width=14, height=10)
  #  GESA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = "  ", main = paste("Primary Sets Assignment - ", doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ") 
  
  dist.matrix <- dist(t(A))
  HC <- hclust(dist.matrix, method="average")
  A <- A[, HC$order]
  A.names <- A.names[HC$order]
  
  dist.matrix <- dist(A)
  HC <- hclust(dist.matrix, method="average")
  A <- A[HC$order,]
  A.row.names <- A.row.names[HC$order]
  
  if (non.interactive.run == FALSE) {
    if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen2, sep=""), sep = "/", collapse="")
      windows(height = height, width = width)
    } else if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  } else {
    if (.Platform$OS.type == "unix") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    } else if (.Platform$OS.type == "windows") {
      filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep=""), sep = "/", collapse="")
      pdf(file=filename, height = height, width = width)
    }
  }
  
  cmap <-  c("#AAAAFF", "#111166")
  
  #   GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  GeneraPipe:::GSEA.HeatMapPlot2(V = t(A), row.names =A.names , col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", sub = paste(doc.string, " - ", phen2, sep=""), xlab=" ", ylab=" ", color.map = cmap) 
  
  text.filename <- paste(directory, paste( doc.string, ".leading.assignment.clustered.", phen2, ".txt", sep=""), sep = "/", collapse="")
  line.list <- c("Gene", A.row.names)
  line.header <- paste(line.list, collapse="\t")
  line.length <- length(A.row.names) + 1
  write(line.header, file = text.filename, ncolumns = line.length)
  write.table(t(A), file=text.filename, append = T, quote=F, col.names= F, row.names=T, sep = "\t")
  
  if (non.interactive.run == FALSE) {  
    if (.Platform$OS.type == "windows") {
      savePlot(filename = filename, type ="jpeg", device = dev.cur())
    } else if (.Platform$OS.type == "unix") {
      dev.off()
    }
  } else {
    dev.off()
  }
  
}