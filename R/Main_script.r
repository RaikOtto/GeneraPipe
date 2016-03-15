#' @export
StartAnalysis = function(
  cel_files_path,
  project,
  msigdb_path,
  output_path = cel_files_path,
  gsea = TRUE,
  kegg_for_heatmap = TRUE
  ){
  
  if (missing(cel_files_path)){
    stop("Need to specify path to input directory.")
  }
  if (missing(project)){
    stop("Need to specify project identifier. Select exisitng project from database or add a new project to database before starting analysis.")
  }
  if (missing(msigdb_path)){
    stop("Need to specify path to MSigDB for GSEA.")
  }
  
  if (project == "GSE29156"){
    cel_files_path = paste(system.file("", package = "GeneraPipe"), "extdata/" , sep = "")
  }
  
  database_path = paste(system.file("", package = "GeneraPipe"), "GeneraPipeDefaultDB.sqlite3", sep = "")
  package_path = system.file("", package = "GeneraPipe")
  
  message( paste( "Running GeneraPipe with project ", project, ".", sep = "" ) )
  #load project parameters
  project_para = data.frame(dplyr::filter(dplyr::tbl(dplyr::src_sqlite(database_path), "projects"), ID == project))
  
  # Set project parameters -----------------------------------------
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
  stat_design  = "contrast"
  heatmap_list_genes_count = 40
  
  # Set default parameters -----------------------------------------
  message("Configurating workspace, loading input data.")  
  paths = GeneraPipe:::set_initial_parameters(
    cel_files_path, 
    output_path,
    database_path,
    project_name,
    cohorts_file,
    kegg_file, 
    package_path,
    msigdb_path)

  pheno_kegg = GeneraPipe:::read_pheno(paths)
  
  # Read in CEL-files -----------------------------------------------
  message("Reading in CEL files." )
  results = GeneraPipe:::read_cel_files(
    cel_files_path,
    chip_type,
    zipped)
  
  # Create cohorts --------------------------------------------------
  message("Creating cohorts for data.")
  results = GeneraPipe:::create_cohorts(
    results,
    chip_type,
    set_ctrl,
    set_case,
    pheno_kegg,
    stat_design)
  
  # Normalization ---------------------------------------------------
  message("Normalizing raw data.")
  results = GeneraPipe:::normalize(
    results,
    chip_type,
    zipped)
  
  # Quality control -------------------------------------------------
  #message("Performing quality control.")
  #GeneraPipe:::quality_control( 
  #  eset, 
  #  raw_data,
  #  phenodata)
  
  # Annotation ------------------------------------------------------
  message("Annotating data.")
  annotation = GeneraPipe:::annotate( 
    results,
    chip_type,
    phenodata)
  
  # Differential expression analysis --------------------------------
  message("Performing differential expression analysis.")
  results = GeneraPipe:::dif_exp_analysis( 
    results,
    chip_type,
    pheno_kegg,
    p_val,
    lfc_exp,
    stat_design)
  
  # Export results --------------------------------------------------
  message(paste("Exporting results to ", paths@results_file_path, ".", sep = ""))
  results = GeneraPipe:::result_preparation(
    paths,
    results,
    annotation,
    chip_type,
    lfc_exp)
  
  # Compute GSEA ----------------------------------------------------
  if (project == "GSE29156"){
    destfile = "~/Output_GeneraPipe/c2.all.v5.1.symbols.gmt"
    download.file("https://raw.githubusercontent.com/janniklas93/GeneraPipeTest/master/c2.all.v5.1.symbols.gmt", destfile, method = "auto")
    paths@msigdb_path = destfile
  }
  if (gsea){
    GeneraPipe:::compute_gsea(
      cel_files_path,
      results,
      set_ctrl,
      set_case,
      annotation,
      paths)
  }
  
  # Create pathway maps ---------------------------------------------
  message(paste("Creating pathway maps, exporting results to ", paths@pathway_maps_path, ".", sep = ""))
  GeneraPipe:::create_pathway_maps(
    results, 
    chip_type, 
    paths,
    annotation,
    pheno_kegg)
  
  # Extract interesting entities ------------------------------------
  message(paste("Computing results for interesting entities, exporting results to ", paths@entities_of_interest_path, ".", sep = "."))
  GeneraPipe:::extract_intr_entities(
    results,
    paths,
    pheno_kegg,
    chip_type,
    annotation)
  
  # Create heatmaps -------------------------------------------------
  if (kegg_for_heatmap){
    message("Creating heatmaps for genes of interest.")
  } else{
    message(paste("Creating heatmaps for the ", heatmap_list_genes_count, " highest differential expressed genes.", sep = ""))
  }
  GeneraPipe:::create_heatmaps(
    results,
    pheno_kegg,
    annotation,
    project_name,
    set_ctrl,
    set_case,
    heatmap_list_genes_count,
    kegg_for_heatmap,
    output_path)
}