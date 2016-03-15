read_pheno = function(paths){
  
  cohorts_file_path = paths@cohorts_file_path
  kegg_file_path    = paths@kegg_file_path
  
  if (GeneraPipe:::strEndsWith(cohorts_file_path, ".csv")){  
    phenodata = read.table( cohorts_file_path , header = TRUE , sep = ",")
  } else {  
    phenodata = read.table( cohorts_file_path , header = TRUE , sep = "\t" )
  }
  phenodata$ID = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", phenodata$ID)
  
  if (GeneraPipe:::strEndsWith(kegg_file_path, ".csv")){  
    keggdata = read.table(kegg_file_path , header = TRUE , sep = ",")
  } else {  
    keggdata = read.table(kegg_file_path , header = TRUE , sep = "\t")
  }
  
  pheno_kegg = list(
    "pheno" = phenodata,
    "kegg"  = keggdata) 
}