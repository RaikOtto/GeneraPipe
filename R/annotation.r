annotate = function(results, chip_type, phenodata){
  
  eset = results@eset

  ### unify all imported ids and check wheather all are necessary in downstream workflow
  
  if ( chip_type == "hgu133plus2" ){  
    hgnc_genes    = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2SYMBOL);
    hgnc_symbols  = hgnc_genes
    hgnc_names    = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GENENAME)
    ensembl_genes = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENSEMBL); 
    entrez_genes  = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENTREZID); 
    uniprot       = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2UNIPROT); 
    go            = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GO); 
    omim          = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2OMIM); 
    enzyme        = BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENZYME); 
      
    hgnc_genes[is.na(hgnc_genes)] = ""
    hgnc_names[is.na(hgnc_names)] = ""
    ensembl_genes[is.na(ensembl_genes)] = ""
    entrez_genes[is.na(entrez_genes)] = ""
    uniprot[is.na(uniprot)] = ""
    go[is.na(go)] = ""
    go[is.na(omim)] = ""
    enzyme[is.na(enzyme)] = ""
  
  } else if ( chip_type == "hgu133a" ){ 
    hgnc_genes    = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aSYMBOL)
    hgnc_symbols  = hgnc_genes
    entrez_genes  = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENTREZID)
    ensembl_genes = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENSEMBL)
    hgnc_names    = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGENENAME)
    uniprot       = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aUNIPROT)
    pathway       = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aPATH)
    go            = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGO)
    omim          = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aOMIM)
    enzyme        = BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENZYME)
  
    entrez_genes[is.na(entrez_genes)] = ""
    hgnc_genes[is.na(hgnc_genes)] = ""
    ensembl_genes[is.na(ensembl_genes)  ] = ""
    hgnc_names[is.na(hgnc_names)  ] = ""
    uniprot[ is.na(uniprot)  ] = ""
    uniprot[ is.na(pathway)  ] = ""
    go[ is.na(go)  ] = ""
    go[ is.na(omim)  ] = ""
    enzyme[ is.na(enzyme)  ] = ""
    
  } else if (chip_type %in% c("pd.hugene.2.0.st")){
    featureData(eset)   = getNetAffx(eset, type = "transcript" )
    hgnc_symbols        = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 2)))
    hgnc_genes          = hgnc_symbols
    hgnc_names          = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 3)))
    ensembl_genes       = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 1)))
  
  } else if (chip_type %in% c("pd.huex.1.0.st.v2")){  
    featureData(eset)   = oligo::getNetAffx(eset, type = "transcript")
    hgnc_symbols        = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 2)))
    hgnc_genes          = hgnc_symbols
    hgnc_names          = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 3)))
    ensembl_genes       = stringr::str_trim(unlist(lapply(featureData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 1)))
    
  } else {
    message("Unknown Chip Type")
    stop()
  }
  
  new("annotation",
    hgnc_genes = hgnc_genes,
    hgnc_symbols = hgnc_symbols,
    hgnc_names = hgnc_names,
    ensembl_genes = ensembl_genes)
}
