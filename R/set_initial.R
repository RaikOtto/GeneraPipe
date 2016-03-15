set_initial_parameters = function(cel_files_path, output_path, database_path, project_name, cohorts_file, kegg_file, package_path, msigdb_path){
  
  # create paths for input and output
  dir.create(paste(output_path, "Output_GeneraPipe", sep = "/"), showWarnings = FALSE)
  cel_files_path              = sub(x = cel_files_path, "/$", "")
  cpdb_file                   = "CPDB_pathways_genes.tab"
  quality_report_path         = paste(output_path, "Output_GeneraPipe/QC_report" , sep = "/" )
  results_file_path           = paste(output_path, "Output_GeneraPipe", paste( "Results", project_name, sep ="_"), sep = "/" )
  name_res_file               = paste(results_file_path, paste( "dif_exp_results", "csv", sep ="."), sep = "/" )
  pathway_maps_path           = paste(results_file_path  , "pathway_maps_dif_exp", sep ="/" )
  gsea_output_path            = paste(paste(output_path, "Output_GeneraPipe/", sep = "/" ), "GSEA_Results", Sys.Date(), sep = "_")
  kegg_file_path              = paste(cel_files_path , kegg_file, sep ="/")
  cpdb_file_path              = paste(package_path, paste( "extdata", cpdb_file, sep ="/" ), sep ="/")
  entities_of_interest_path   = paste(results_file_path, "Entities_of_interest", sep ="/")
  genes_of_interest_file_path = paste(entities_of_interest_path, "genes_of_interest.xls", sep ="/") 
  user_folder                 = as.character( system("echo $HOME", intern = TRUE) )
  cohorts_file_path           = paste(cel_files_path, cohorts_file, sep ="/" )
  
  new("paths",
      quality_report_path         = quality_report_path,
      results_file_path           = results_file_path,
      name_res_file               = name_res_file,
      pathway_maps_path           = pathway_maps_path,
      gsea_output_path            = gsea_output_path,
      kegg_file_path              = kegg_file_path,
      cpdb_file_path              = cpdb_file_path,
      entities_of_interest_path   = entities_of_interest_path,
      genes_of_interest_file_path = genes_of_interest_file_path,
      user_folder                 = user_folder,
      cohorts_file_path           = cohorts_file_path,
      package_path                = package_path,
      msigdb_path                 = msigdb_path
  )
  
  # create paths for input and output
  #dir.create(paste(output_path, "Output_GeneraPipe", sep = "/"), showWarnings = FALSE)
  #cel_files_path                  = sub(x = cel_files_path, "/$", "")
  #cpdb_file                       = "CPDB_pathways_genes.tab"
  #env$quality_report_path         = paste(output_path, "Output_GeneraPipe/QC_report" , sep = "/" )
  #env$results_file_path           = paste(output_path, "Output_GeneraPipe", paste( "Results", project_name, sep ="_"), sep = "/" )
  #env$name_res_file               = paste(env$results_file_path, paste( "dif_exp_results", "csv", sep ="."), sep = "/" )
  #env$pathway_maps_path           = paste(env$results_file_path  , "pathway_maps_dif_exp", sep ="/" )
  #env$gsea_output_path            = paste(output_path, "Output_GeneraPipe/GSEA_Results", sep = "/" )
  #env$kegg_file_path              = paste(cel_files_path , kegg_file, sep ="/")
  #env$cpdb_file_path              = paste(package_path, paste( "extdata", cpdb_file, sep ="/" ), sep ="/")
  #env$entities_of_interest_path   = paste(env$results_file_path, "Entities_of_interest", sep ="/")
  #env$genes_of_interest_file_path = paste(env$entities_of_interest_path, "genes_of_interest.xls", sep ="/") 
  #env$user_folder                 = as.character( system("echo $HOME", intern = TRUE) )
  #cohorts_file_path               = paste(cel_files_path, cohorts_file, sep ="/" )
  
  #return(phenodata)
}