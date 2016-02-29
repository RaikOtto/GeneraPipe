dif_exp_analysis = function( eset, chip_type, phenodata, p_val, lfc_exp ){
  
  if (chip_type %in% c("hgu133plus2", "hgu133a")){
    eset = eset[! gdata::startsWith(rownames(eset), "AFFX-"),]
  }

  eset = eset[, which(colnames(eset) %in% phenodata$ID)]
  index = which(colnames(phenodata) == env$cohorts_type )
  pData(eset)$Group = phenodata[match(colnames(eset), phenodata$ID, nomatch = 0), index]
  eDatSet = eset

  if (env$var_filter){

    exprs(eDatSet) = exprs(genefilter::varFilter(eDatSet) )

  }

  #source("Src/cohort_creation.r")
  fit = limma::lmFit(eDatSet[, ], env$design)

  if (env$stat_design == "contrast"){

      #fit = lmFit( eDatSet[ ,c( index_case, index_ctrl )  ], design )
      #fit = lmFit( eDatSet[ ,index_cohorts_vec  ], design )
      cont.matrix = limma::makeContrasts(contrast = CASE - CTRL, levels = env$design)
      fit = limma::contrasts.fit(fit, cont.matrix)
      fit = limma::eBayes(fit)
      #volc_all = topTable(fit, adjust="fdr", sort.by="B", number = 50000)
      volc_all = limma::topTable(fit, coef = "contrast", number  = nrow(eDatSet), adjust  ="none", p.value = 1, lfc = 0)

  } else {

      fit = limma::eBayes(fit)
      volc_all = limma::topTable(fit, number  = nrow(eDatSet), adjust  ="none", p.value = 1, lfc = 0)
  }

  dir.create(env$results_file_path, showWarnings = F)

  topall = limma::topTable(fit, coef = "contrast", number  = nrow( eDatSet ), adjust  = "none", p.value = p_val, lfc = lfc_exp)

  if ((dim(topall)[1] == 0) & (dim(topall)[2] == 0)){
    stop("Topall has dimension zero")

  }

  topall$logFC = round(topall$logFC,2)
  topall$AveExpr = round(topall$AveExpr,2)
  topall$t = round( topall$t, 2 )
  topall$B = round( topall$B, 2 )

  topall = topall[ abs(topall$logFC) >= lfc_exp  ,]

  message( c( "Amount probes higher in Case cohort: ", sum( topall$logFC >= lfc_exp ) ) )
  message( c( "Amount probes lower in Case cohort: " , sum( topall$logFC < lfc_exp ) ) )
  
  env$eset = eset
  
  env$volc_all= volc_all
  return(topall)
}