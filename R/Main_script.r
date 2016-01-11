#' @export
start_analysis = function(cel_files_path = "", project = "GSE14407"){
  
  package_path = system.file("", package = "GeneraPipe")
  database_path = paste(package_path, "GeneraPipeDefaultDB.sqlite3", sep = "/")
  
  #load project parameters
  project_para = data.frame(dplyr::filter(tbl(src_sqlite(database_path), "projects"), ID == project))
  
  #extract parameters
  chip_type = project_para$chip_type
  if (project_para$zipped == 1){
    zipped = TRUE
  } else{
    zipped = FALSE
  }
  
  parse_cel_files(cel_files_path, chip_type, zipped)

}
