read_cel_files = function(cel_files_path, chip_type, zipped){
  
  if (chip_type %in% c("hgu133plus2", "hgu133a")){
    celFiles = oligoClasses::list.celfiles(cel_files_path, full = TRUE, listGzipped = zipped)
    raw_data = affy::read.affybatch( filenames = celFiles)
    
    #results = new("results",
        #raw_data1 = raw_data)
  
  } else if (chip_type %in% c("pd.huex.1.0.st.v2", "pd.hugene.2.0.st")){
    celFiles = oligoClasses::list.celfiles( cel_files_path, full = TRUE, listGzipped = zipped)
    raw_data = oligo::read.celfiles(filenames = celFiles)
    
    #results = new("results",
        #raw_data2 = raw_data)
  
  } else {
    message(c("Unknown Chip Type: ", chip_type))
    stop()
  }
  
  #results
  new("results",
      raw_data = raw_data,
      eset = ExpressionSet(),
      topall = data.frame(),
      volc_all = data.frame(),
      cohorts_vec = NA_character_
  )
}
