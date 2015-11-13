### pathway mapping

library("stringr")

kegg_t = read.table( kegg_file_path, header =T , sep ="\t" )
cpdb_t = read.table( cpdb_file_path, header =T , sep ="\t", fill = T )
cpdb_ident = str_replace(cpdb_t$external_id, "path:", "")

interesting_pathways_mapping = match( str_trim( kegg_t$HSA_ID ), str_trim( cpdb_ident ) )
interesting_pathways_table = cpdb_t[ interesting_pathways_mapping ,]

if ( chip_type == "hgu133plus2" ){
  
  if ( multi_probe ){
    
    ### yet to be done
    message("Multiprobe annotation not implemented yet, System aborting")
    quit()
    
  } else {
    
    hgnc_genes = mget( rownames(eset), hgu133plus2SYMBOL ); hgnc_genes[ is.na(hgnc_genes)  ] = ""
    ensembl_genes = mget( rownames(eset), hgu133plus2ENSEMBL ); ensembl_genes[ is.na(ensembl_genes)  ] = ""
    entrez_genes = mget( rownames(eset), hgu133plus2ENTREZID ); entrez_genes[ is.na(entrez_genes)  ] = ""
    uniprot = mget( rownames(eset), hgu133plus2UNIPROT ); uniprot[ is.na(uniprot)  ] = ""
    #pfam   = select( rownames(eset), hgu133plus2PFAM )
    go = mget( rownames(eset), hgu133plus2GO ); go[ is.na(go)  ] = ""
    omim = mget( rownames(eset), hgu133plus2OMIM ); go[ is.na(omim)  ] = ""
    enzyme = mget( rownames(eset), hgu133plus2ENZYME ); enzyme[ is.na(enzyme)  ] = ""
    
  }
  
} else if ( chip_type == "hgu133a" ){
  
  if ( multi_probe ){
    
    mapWithMultiProbes_entrez = toggleProbes( hgu133aENTREZID, "all")
    mapWithMultiProbes_symbols = toggleProbes( hgu133aSYMBOL, "all")
    
    hgnc_genes = mget( rownames( eset ), mapWithMultiProbes_symbols )
    entrez_genes = mget( rownames( eset ), mapWithMultiProbes_entrez )
    
  } else {
    
    hgnc_genes = mget( rownames(eset), hgu133aSYMBOL ); hgnc_genes[ is.na(hgnc_genes)  ] = ""
    entrez_genes = mget( rownames(eset), hgu133aENTREZID ); entrez_genes[ is.na(entrez_genes)  ] = ""
  }  
  
  ensembl_genes = mget( rownames(eset), hgu133aENSEMBL ); ensembl_genes[ is.na(ensembl_genes)  ] = ""
  hgnc_names = mget( rownames(eset), hgu133aGENENAME ); hgnc_genes[ is.na(hgnc_names)  ] = ""
  uniprot = mget( rownames(eset), hgu133aUNIPROT ); uniprot[ is.na(uniprot)  ] = ""
  pathway = mget( rownames(eset), hgu133aPATH ); uniprot[ is.na(pathway)  ] = ""
  #pfam   = select( rownames(eset), PFAM )
  go = mget( rownames(eset), hgu133aGO ); go[ is.na(go)  ] = ""
  omim = mget( rownames(eset), hgu133aOMIM ); go[ is.na(omim)  ] = ""
  enzyme = mget( rownames(eset), hgu133aENZYME ); enzyme[ is.na(enzyme)  ] = ""
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st") ){
  
  featureData( eset  ) = getNetAffx( eset, type = "transcript" )
  split_fun = function( entry, pos ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[pos] ) } else{ return( "" ) } }
  hgnc_symbols = str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN=split_fun, 2 ) ) )
  hgnc_genes = hgnc_symbols
  hgnc_names = str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN=split_fun, 3 ) ) )
  ensembl_genes = str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN=split_fun, 1 ) ) )
  
} else if ( chip_type %in% c( "pd.huex.1.0.st.v2" ) ){  
  
  featureData( eset  ) = getNetAffx( eset, type = "transcript" )
  split_fun = function( entry, pos ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[pos] ) } else{ return( "" ) } }
  hgnc_symbols = str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN=split_fun, 2 ) ) )
  
} else {
  
  message("Unknown Chip Type")
  stop()
}

if (integrate_drug_data){
  
  
  index_drug = match( drug_type, colnames( phenodata ), nomatch = 0 )
  
  if (index_drug == 0){
    
    message("Could not find drug in cohorts file")
    quit()
  }
  
  mapping = match( names(cohorts_vec), phenodata$ID, nomatch = 0)
  drug_data = phenodata[ mapping  ,index_drug ]
  names(drug_data) = names(cohorts_vec)
}
