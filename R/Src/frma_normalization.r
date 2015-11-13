print("Running step 4b: fRMA Normalizing"  )

setwd(frma_path)

if ( chip_type %in% c( "hgu133a", "hgu133plus2" ) ){
  
  eset = threestep( raw_data, background.method = "RMA.2", normalize.method="quantile", summary.method="median.polish")
  
} else if ( chip_type %in% c( "pd.huex.1.0.st.v2" ) ){

  eset = rma( raw_data, target='extended', normalize = T)
  print("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st" ) ){
  
  eset = rma( raw_data, target='core', normalize = T)
  print("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else {
  
  print(c("Unknown Chip Type: ",chip_type))
  stop()
}

#pData(eset) = phenodata
if (! zipped ){
  mapping_cohort_p = match( rownames(pData(eset)), names(cohorts_vec) , nomatch = 0 )
  mapping_cohort_c = match( names(cohorts_vec), rownames(pData(eset)) , nomatch = 0 )
  
} else{
  
  mapping_cohort_p = match( rownames(pData(eset)), paste0( names(cohorts_vec), ".gz" ), nomatch = 0 )
  mapping_cohort_c = match( paste0( names(cohorts_vec), ".gz"), rownames(pData(eset)) , nomatch = 0 )
  
}

pData(eset)$Cohort = cohorts_vec[ mapping_cohort_p ]
