read_cel_files = function( cel_files_path, chip_type, zipped ){
  
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
  
  } else if ( chip_type %in% c( "pd.huex.1.0.st.v2", "pd.hugene.2.0.st" ) ){
  
    celFiles = oligoClasses::list.celfiles( cel_files_path, full = TRUE, listGzipped = zipped)
    raw_data = oligo::read.celfiles( filenames = celFiles )
  
  } else {
  
    message( c("Unknown Chip Type: ", chip_type) )
    stop()
  }
  
  return(raw_data)
}
