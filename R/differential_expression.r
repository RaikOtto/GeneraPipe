dif_exp_analysis = function(results, chip_type, pheno_kegg, p_val, lfc_exp, stat_design){
  
  eset = results@eset
  design = results@design
  phenodata = pheno_kegg$pheno
  
  eset = eset[, which(colnames(eset) %in% phenodata$ID)]
  index = which(colnames(phenodata) == "Group")
  pData(eset)$Group = phenodata[match(colnames(eset), phenodata$ID, nomatch = 0), index]

  ## if var_filter stays in pipeline, muss eset für diese analyse einer anderen variable zugewiesen werden (eDatSet = eset))
  #if (env$var_filter){

   # exprs(eset) = exprs(genefilter::varFilter(eset) )

  #}

  fit = limma::lmFit(eset[, ], design)

  if (stat_design == "contrast"){
      cont.matrix = limma::makeContrasts(contrast = CASE - CTRL, levels = design)
      fit = limma::contrasts.fit(fit, cont.matrix)
      fit = limma::eBayes(fit)
      volc_all = limma::topTable(fit, coef = "contrast", number  = nrow(eset), adjust  = "none", p.value = 1, lfc = 0)

  } else {
      fit = limma::eBayes(fit)
      volc_all = limma::topTable(fit, number  = nrow(eset), adjust  ="none", p.value = 1, lfc = 0)
  }

  topall = limma::topTable(fit, coef = "contrast", number  = nrow(eset), adjust  = "none", p.value = p_val, lfc = lfc_exp)
  index_topall = match(row.names(topall), row.names(eset))

  if ((dim(topall)[1] == 0) & (dim(topall)[2] == 0)){
    stop("Topall has dimension zero")
  }

  topall$logFC   = round(topall$logFC,2)
  topall$AveExpr = round(topall$AveExpr,2)
  topall$t       = round(topall$t, 2)
  topall$B       = round(topall$B, 2)
  topall         = topall[abs(topall$logFC) >= lfc_exp, ]

  message(c("Amount probes higher in Case cohort: ", sum(topall$logFC >= lfc_exp)))
  message(c("Amount probes lower in Case cohort: ", sum(topall$logFC < lfc_exp)))
  
  results@eset         = eset
  results@volc_all     = volc_all
  results@topall       = topall
  results@index_topall = index_topall
  results
}