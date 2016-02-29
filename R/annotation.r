annotate = function( eset, chip_type, phenodata ){

  if ( chip_type == "hgu133plus2" ){
  
    if ( env$multi_probe ){   
      ### yet to be done
      message("Multiprobe annotation not implemented yet, System aborting")
      quit()
    } else {  
      env$hgnc_genes    = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2SYMBOL); 
      env$ensembl_genes = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENSEMBL); 
      env$entrez_genes  = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENTREZID); 
      env$uniprot       = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2UNIPROT); 
      env$go            = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GO); 
      env$omim          = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2OMIM); 
      env$enzyme        = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENZYME); 
      
      env$hgnc_genes[ is.na(env$hgnc_genes)  ] = ""
      env$ensembl_genes[ is.na(env$ensembl_genes)  ] = ""
      env$entrez_genes[ is.na(env$entrez_genes)  ] = ""
      env$uniprot[ is.na(env$uniprot)  ] = ""
      env$go[ is.na(env$go)  ] = ""
      env$go[ is.na(env$omim)  ] = ""
      env$enzyme[ is.na(env$enzyme)  ] = ""
    }
  
  } else if ( chip_type == "hgu133a" ){
  
    if ( multi_probe ){
    
      mapWithMultiProbes_entrez  = AnnotationDbi::toggleProbes(hgu133a.db::hgu133aENTREZID, "all")
      mapWithMultiProbes_symbols = AnnotationDbi::toggleProbes(hgu133a.db::hgu133aSYMBOL, "all")
    
      env$hgnc_genes   = BiocGenerics::mget(rownames(eset), mapWithMultiProbes_symbols)
      env$entrez_genes = BiocGenerics::mget(rownames(eset), mapWithMultiProbes_entrez)
    
    } else {
    
      env$hgnc_genes   = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aSYMBOL )
      env$entrez_genes = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENTREZID )
      
      env$entrez_genes[is.na(entrez_genes)] = ""
      env$hgnc_genes[is.na(hgnc_genes)] = ""
    }  
  
    env$ensembl_genes = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENSEMBL)
    env$hgnc_names    = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGENENAME)
    env$uniprot       = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aUNIPROT)
    env$pathway       = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aPATH)
    env$go            = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGO)
    env$omim          = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aOMIM)
    env$enzyme        = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENZYME)
  
    env$ensembl_genes[ is.na(env$ensembl_genes)  ] = ""
    env$hgnc_genes[ is.na(env$hgnc_names)  ] = ""
    env$uniprot[ is.na(env$uniprot)  ] = ""
    env$uniprot[ is.na(env$pathway)  ] = ""
    env$go[ is.na(env$go)  ] = ""
    env$go[ is.na(env$omim)  ] = ""
    env$enzyme[ is.na(env$enzyme)  ] = ""
    
  } else if ( chip_type %in% c( "pd.hugene.2.0.st") ){
  
    featureData(eset) = getNetAffx(eset, type = "transcript" )
    env$hgnc_symbols        = stringr::str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN = GeneraPipe::split_fun, 2 ) ) )
    env$hgnc_genes          = env$hgnc_symbols
    env$hgnc_names          = stringr::str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN = GeneraPipe::split_fun, 3 ) ) )
    env$ensembl_genes       = stringr::str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN = GeneraPipe::split_fun, 1 ) ) )
  
  } else if ( chip_type %in% c( "pd.huex.1.0.st.v2" ) ){  
  
    featureData(eset)      = oligo::getNetAffx(eset, type = "transcript")
    env$hgnc_symbols = stringr::str_trim( unlist( lapply( featureData( eset )$geneassignment, FUN = GeneraPipe::split_fun, 2 ) ) )
    env$hgnc_genes          = env$hgnc_symbols
    
  } else {
  
    message("Unknown Chip Type")
    stop()
  }

  if (env$integrate_drug_data){
    index_drug = match(drug_type, colnames(phenodata), nomatch = 0)
    if (index_drug == 0){
    
      message("Could not find drug in cohorts file")
      quit()
    }
    mapping = match(names(env$cohorts_vec), phenodata$ID, nomatch = 0)
    drug_data = phenodata[mapping  ,index_drug]
    names(drug_data) = names(env$cohorts_vec)
  }
  
  return(eset)
}
