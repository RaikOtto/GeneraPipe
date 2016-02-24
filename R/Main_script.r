#' @export
start_analysis = function(cel_files_path = "", output_path = cel_files_path, project = "GSE14407"){
  
  # cel_files_path = "~/Generic_mRNA_Expression_Pipeline/Project_files/GSE14407_RAW/Input/"
  # output_path = "~/Documents"
  package_path = system.file("", package = "GeneraPipe")
  database_path = "~/GeneraPipe/inst/GeneraPipeDefaultDB.sqlite3"
  #database_path = paste(package_path, "GeneraPipeDefaultDB.sqlite3", sep = "/")
  
  message( paste( "Running GeneraPipe with project ", project, ".", sep = "" ) )
  #load project parameters
  project_para = data.frame(dplyr::filter(tbl(src_sqlite(database_path), "projects"), ID == project))
  
  #set project parameters
  message( "Setting project parameters." )
  package_env$project_name = project_para$project_name
  chip_type                = project_para$chip_type
  package_env$kegg_file    = project_para$kegg_file
  package_env$cohorts_file = project_para$cohorts_file
  set_case                 = project_para$set_case
  set_ctrl                 = project_para$set_ctrl
  package_env$p_val        = project_para$p_val
  package_env$lfc_exp      = project_para$lfc_exp
  zipped                   = as.logical(project_para$zipped)
  
  # set default parameters
  message( "Configurating workspace, loading input data.")  
  GeneraPipe::set_initial_parameters(cel_files_path, output_path, database_path)
  
  # read in CEL files
  message( "Reading in CEL files." )
  raw_data = GeneraPipe::read_cel_files(cel_files_path, chip_type, project_name, zipped)
  
  # create cohorts
  message( "Creating cohorts for data.")
  raw_data = GeneraPipe::create_cohorts( raw_data, chip_type, set_ctrl, set_case)
  
  # Normalization
  message("Normalizing raw data.")
  eset = GeneraPipe::normalize(raw_data = raw_data, chip_type = chip_type, zipped = zipped)
  
  # Annotation
  message( "Annotating data." )
  GeneraPipe::annotate( eset, chip_type )
  
}
