#suppressMessages(library("arrayQualityMetrics"))
quality_control = function( eset, raw_data, phenodata ){

  if (! (quality_control_only)){
  
    mapping = match(colnames(eset), phenodata$ID)
  
    for (elem in colnames(phenodata)){
    
      if ((elem != "ID") && (! (elem %in% colnames(pData(eset))))){ 
      
        index_elem = which(colnames(phenodata) == elem)
        pData(eset)$elem = phenodata[mapping, index_elem]
        colnames(pData(eset))[length(colnames(pData(eset)))] = elem
      }
    }
    pData(eset)$Cohort = cohorts_vec[match(rownames(pData(eset)), names(env$cohorts_vec) ) ]
  
  } else {
  
    mapping_eset_pheno = match(as.character(colnames(eset)), as.character(phenodata$ID))
    pData(eset)$Group = phenodata$Group[mapping_eset_pheno]
  }

  pData(raw_data) = pData(eset)

  message( "Running Quality Metrics"  )
  arrayQualityMetrics( raw_data, intgroup = "Cohort", outdir = env$quality_report_path, force = TRUE, showWarnings = FALSE)
  #arrayQualityMetrics( raw_data, force = T)
}