annotate = function(results, chip_type, phenodata){
  
  eset = results@eset
  
  #if (chip_type %in% c("hgu133plus2", "hgu133a", "drosophila2")){
  #  eset = eset[! gdata::startsWith(rownames(eset), "AFFX-"),]
  #}
  
  ### unify all imported ids and check wheather all are necessary in downstream workflow
  
  if (chip_type == "hgu133plus2"){  
    hgnc_genes    = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2SYMBOL))
    hgnc_symbols  = hgnc_genes
    hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GENENAME))
    ensembl_genes = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENSEMBL))
    entrez_genes  = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENTREZID))
    uniprot       = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2UNIPROT))
    go            = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2GO))
    omim          = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2OMIM))
    enzyme        = unlist(BiocGenerics::mget(rownames(eset), hgu133plus2.db::hgu133plus2ENZYME))
      
    hgnc_genes[is.na(hgnc_genes)] = ""
    hgnc_names[is.na(hgnc_names)] = ""
    ensembl_genes[is.na(ensembl_genes)] = ""
    entrez_genes[is.na(entrez_genes)] = ""
    uniprot[is.na(uniprot)] = ""
    go[is.na(go)] = ""
    go[is.na(omim)] = ""
    enzyme[is.na(enzyme)] = ""
  
  } else if (chip_type == "hgu133a"){ 
    hgnc_genes    = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aSYMBOL))
    hgnc_symbols  = hgnc_genes
    entrez_genes  = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENTREZID))
    ensembl_genes = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENSEMBL))
    hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGENENAME))
    uniprot       = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aUNIPROT))
    pathway       = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aPATH))
    go            = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aGO))
    omim          = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aOMIM))
    enzyme        = unlist(BiocGenerics::mget(rownames(eset), hgu133a.db::hgu133aENZYME))
  
    entrez_genes[is.na(entrez_genes)] = ""
    hgnc_genes[is.na(hgnc_genes)] = ""
    ensembl_genes[is.na(ensembl_genes)  ] = ""
    hgnc_names[is.na(hgnc_names)  ] = ""
    uniprot[ is.na(uniprot)  ] = ""
    uniprot[ is.na(pathway)  ] = ""
    go[ is.na(go)  ] = ""
    go[ is.na(omim)  ] = ""
    enzyme[ is.na(enzyme)  ] = ""
    
  } else if (chip_type %in% c("drosophila2")){
    
    hgnc_genes    = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2SYMBOL))
    hgnc_symbols  = hgnc_genes
    entrez_genes  = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2ENTREZID))
    ensembl_genes = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2ENSEMBL))
    hgnc_names    = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2GENENAME))
    uniprot       = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2UNIPROT))
    pathway       = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2PATH))
    go            = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2GO))
    enzyme        = unlist(BiocGenerics::mget(rownames(eset), drosophila2.db::drosophila2ENZYME))
 
    hgnc_genes[is.na(hgnc_genes)] = ""
    entrez_genes[is.na(entrez_genes)] = ""
    ensembl_genes[is.na(ensembl_genes)] = ""
    hgnc_names[is.na(hgnc_names)] = ""
    uniprot[is.na(uniprot)] = ""
    uniprot[is.na(pathway)] = ""
    go[is.na(go)] = ""
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
