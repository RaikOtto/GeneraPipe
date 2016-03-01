#' @export
start_analysis = function(
  cel_files_path = paste(system.file("", package = "GeneraPipe"), "extdata", sep = "/"),
  output_path = "",
  database_path = paste(system.file("", package = "GeneraPipe"), "GeneraPipeDefaultDB.sqlite3", sep = "/"),
  project = "GSE29156",
  gene_set_database_path = ""
  ){
  
  #cel_files_path = "~/Documents/GeneraPipe_GSE29156/Input"
  #output_path = "~/Documents"
  package_path = system.file("", package = "GeneraPipe")
  #database_path = "~/GeneraPipe/inst/GeneraPipeDefaultDB.sqlite3"
  #database_path = paste(package_path, "GeneraPipeDefaultDB.sqlite3", sep = "/")
  
  message( paste( "Running GeneraPipe with project ", project, ".", sep = "" ) )
  #load project parameters
  project_para = data.frame(dplyr::filter(dplyr::tbl(dplyr::src_sqlite(database_path), "projects"), ID == project))
  
  #set project parameters
  message("Setting project parameters.")
  project_name = project_para$project_name
  chip_type    = project_para$chip_type
  kegg_file    = project_para$kegg_file
  cohorts_file = project_para$cohorts_file
  set_case     = project_para$set_case
  set_ctrl     = project_para$set_ctrl
  p_val        = project_para$p_val
  lfc_exp      = project_para$lfc_exp
  zipped       = as.logical(project_para$zipped)
  
  # set default parameters
  message( "Configurating workspace, loading input data.")  
  phenodata = GeneraPipe:::set_initial_parameters(cel_files_path, output_path, database_path, project_name, cohorts_file, kegg_file, package_path)
  
  # read in CEL files
  message( "Reading in CEL files." )
  raw_data = GeneraPipe:::read_cel_files(cel_files_path, chip_type, zipped)
  
  # create cohorts
  message( "Creating cohorts for data.")
  raw_data = GeneraPipe:::create_cohorts( raw_data, chip_type, set_ctrl, set_case, phenodata)
  phenodata = env$phenodata
  
  # normalization
  message("Normalizing raw data.")
  eset = GeneraPipe:::normalize(raw_data = raw_data, chip_type = chip_type, zipped = zipped)
  
  # annotation
  message("Annotating data.")
  eset = GeneraPipe:::annotate( eset, chip_type, annotate )
  
  # differential expression analysis
  message("Performing differential expression analysis.")
  topall = GeneraPipe:::dif_exp_analysis( eset, chip_type, phenodata, p_val, lfc_exp )
  eset = env$eset
  volc_all = env$volc_all
  
  # export results
  message(paste("Exporting results to ", env$results_file_path, ".", sep = ""))
  topall_res = GeneraPipe:::result_preparation( eset, topall, chip_type, lfc_exp )
  eset = env$eset
  
  # compute GSEA
  #if (env$use_gsea){
  #  GeneraPipe:::compute_gsea( cel_files_path, eset, set_case, set_ctrl, gene_set_database_path )
  #}
  
  # create pathway maps
  message(paste("Creating pathway maps, exporting results to ", env$pathway_maps_path, ".", sep = ""))
  GeneraPipe:::create_pathway_maps( topall_res, eset, volc_all, chip_type, package_path )
  
  # extract interesting entities
  message(paste("Computing results for interesting entities, exporting results to ", env$entities_of_interest_path, ".", sep = "."))
  GeneraPipe:::extract_intr_entities( eset, chip_type )
  
  # create heatmaps
  if (env$create_heatmaps_genes_of_interest){
    if (env$use_kegg_for_heatmap){
      message("Creating heatmaps for genes of interest.")
    } else{
      message(paste("Creating heatmaps for the ", env$heatmap_list_genes_count, " highest differential expressed genes.", sep = ""))
    }
    GeneraPipe:::create_heatmaps( eset, topall_res, project_name, set_ctrl, set_case )
  }
}