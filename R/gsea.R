compute_gsea = function(cel_files_path, results, set_ctrl, set_case, annotation, paths){
  
  gsea_output_path = paths@gsea_output_path
  msigdb_path = paths@msigdb_path
  eset       = results@eset
  index_case = results@index_case
  index_ctrl = results@index_ctrl
  hgnc_symbols = annotation@hgnc_symbols
  
  if (! file.exists(paste(cel_files_path, "ExpressionSet.gct", sep = "/"))){
    message("Creating ExpressionSet.gct file for GSEA.")
    probe_ids = hgnc_symbols
    descr = rep("NA", length( probe_ids))
    data = as.data.frame(exprs(eset))
    expression_out = data.frame("NAMES" = probe_ids, "DESCRIPTION" = descr, data)
  
    index_case_gsea = sapply(index_case, function(x) x + 2)
    index_ctrl_gsea = sapply(index_ctrl, function(x) x + 2)
  
    expression_out = expression_out[expression_out$NAMES != "", ]
    names_dupli = expression_out$NAMES[which( duplicated( expression_out$NAMES))]
    names_dupli = unique(names_dupli)
  
    reduced = list()
  
    if( length(names_dupli) > 0){
    
      # collapsing duplicate gene symbols by logFC
      index_dupli = which(expression_out$NAMES %in% names_dupli)
      for (i in 1:length(names_dupli)){
      
        index = which(expression_out$NAMES == names_dupli[i])
        current = expression_out[index, ]
      
        exprs_case_gsea = rowMeans(current[, index_case_gsea])
        exprs_ctrl_gsea = rowMeans(current[, index_ctrl_gsea])
        dif_exp_gsea    = exprs_case_gsea - exprs_ctrl_gsea
      
        index_max = which.max(abs(dif_exp_gsea))
        temp      = as.data.frame(current[index_max, ])
        reduced   = c(reduced, list(temp))
      }
    
      reduced = do.call("rbind", reduced)
      expression_out = expression_out[-index_dupli, ]
      expression_out = rbind(expression_out, reduced)
    }
  
    # write ExpressionSet.gct file in Input directory
    fileConn = file(paste(cel_files_path, "ExpressionSet.gct", sep = "/" ))
    writeLines(c("#1.2", paste(length(expression_out$NAMES), length(expression_out) - 2, sep = "\t")), fileConn)
    close(fileConn)
    write.table(expression_out, file = paste(cel_files_path, "ExpressionSet.gct", sep ="/"), sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  
  } else{
    message(paste("There already exists ExpressionSet.gct file in ", cel_files_path, ".", " Using this for GSEA.", sep = ""))
  }

  if (! file.exists(paste(cel_files_path, "phenotypes_GSEA.cls", sep = "/"))){
    message( "Creating phenotypes_GSEA.cls file for GSEA." )
    if(! exists("expression_out")){
      probe_ids = hgnc_symbols
      descr = rep("NA", length( probe_ids))
      data = as.data.frame(exprs(eset))
      expression_out = data.frame("NAMES" = probe_ids, "DESCRIPTION" = descr, data)
    }
    fileConn = file(paste(cel_files_path, "phenotypes_GSEA.cls", sep = "/"))
    phenoLabels = rep(0, length(expression_out) - 2)
    case_GSEA = unlist(strsplit(set_case, " "))
    case_GSEA = case_GSEA[1]
    ctrl_GSEA = unlist(strsplit(set_ctrl, " "))
    ctrl_GSEA = ctrl_GSEA[1]
    phenoLabels[index_case] = case_GSEA
    phenoLabels[index_ctrl] = ctrl_GSEA
    phenoLabels = paste(phenoLabels, collapse = " ")
    if (phenoLabels[1] == set_case){
      pheno_order = paste("#", case_GSEA, ctrl_GSEA, sep = " ")
    } else {
      pheno_order = paste("#", ctrl_GSEA, case_GSEA, sep = " ")
    }
    writeLines(c(paste(length(expression_out) - 2 , "2", "1", sep = " "), pheno_order, phenoLabels), fileConn)
    close(fileConn)
    #write.table(phenoLabels, file = paste( cel_files_path, "phenotypes_GSEA.cls", sep ="/" ), sep = " ", col.names = FALSE, row.names = F, append=TRUE)
  } else{
    message(paste("There already exists phenotypes_GSEA.cls file in ", cel_files_path, ".", " Using this for GSEA.", sep = ""))
  }

  unlink(paste(gsea_output_path, "*", sep = "/"))

  GeneraPipe:::GSEA(                                                         # Input/Output Files :-------------------------------------------
                                                                           input.ds              = paste(cel_files_path, "ExpressionSet.gct", sep ="/"),               # Input gene expression Affy dataset file in RES or GCT format
                                                                           input.cls             = paste(cel_files_path, "phenotypes_GSEA.cls", sep ="/"),               # Input class vector (phenotype) file in CLS format
                                                                           gs.db                 = msigdb_path,           # Gene set database in GMT format
                                                                           output.directory      = gsea_output_path,            # Directory where to store output and results (default: "")
                                                                           #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
                                                                           doc.string            = "Benign_Stroma_vs_Malignant_Stroma",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                                                                           non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                                                                           reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                                                                           nperm                 = 1000,            # Number of random permutations (default: 1000)
                                                                           weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                                                                           nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                                                                           fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                                                                           fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                                                                           topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                                                                           adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                                                                           gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                                                                           gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                                                                           reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                                                                           preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                                                                           random.seed           = 111,             # Random number generator seed. (default: 123456)
                                                                           perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                                                                           fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                                                                           replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                                                                           save.intermediate.results = TRUE,        # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                                                                           OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                                                                           use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
  )
  #--------------------------------------------------------------------------------------------------------------------------------------------------

  # Overlap and leading gene subset assignment analysis of the GSEA results                                                                           
  GeneraPipe:::GSEA.Analyze.Sets(
    directory = gsea_output_path,        # Directory where to store output and results (default: "")
    topgs     = 20,                     # number of top scoring gene sets used for analysis
    height    = 16,
    width     = 16
  )
  
}