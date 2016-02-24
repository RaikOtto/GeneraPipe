read_cel_files = function(cel_files_path, chip_type, project_name, zipped){
  
  if ( chip_type == "HumanHT-12.v4" ){
  
    raw_data = getGEO(GEO_Id, GSEMatrix = TRUE)
    if ( length( raw_data ) > 1 )
      idx = grep("GPL10558", attr(raw_data, "names"))
    else
      idx = 1
    raw_data = raw_data[[idx]]
  
  } else if ( chip_type %in% c( "hgu133plus2", "hgu133a" ) ){
  
    if ( zipped ){
      celFiles = affy::list.celfiles( cel_files_path, full = TRUE)
    } else {
      celFiles = affy::list.celfiles( cel_files_path, full = TRUE)
    }
  
    raw_data = affy::read.affybatch( filenames = celFiles )
    #excluded_files = c()
  
    if( package_env$project_name == "Plus2_Runx"  ){
      excluded_files = c( 1,2,3,21,94,115,116,117,118,119,120,121,122,123,124,125,129,138,153,154,155,156,157,158,159,160,161,164,165,174,175,176,179)
      raw_data = raw_data[ -( excluded_files )  ] # exclude bad samples
      phenodata = phenodata[ match( phenodata$ID, colnames(raw_data), nomatch = 0 ),]
    }
  
  } else if ( chip_type %in% c( "pd.huex.1.0.st.v2", "pd.hugene.2.0.st" ) ){
  
    celFiles = oligoClasses::list.celfiles( cel_files_path, full = TRUE, listGzipped = zipped)
    raw_data = oligo::read.celfiles( filenames = celFiles )
  
  } else {
  
    message( c("Unknown Chip Type: ", chip_type) )
    stop()
  }
  
  return(raw_data)
}
