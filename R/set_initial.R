env <- new.env(parent = emptyenv())

set_initial_parameters = function( cel_files_path, output_path, database_path, project_name, cohorts_file, kegg_file, package_path ){
  
  # set default parameters, extracted from db
  default_pm = data.frame(dplyr::filter( dplyr::tbl(dplyr::src_sqlite(database_path), "default_pm")))
  env$use_gsea                          = as.logical(default_pm$use_gsea)
  env$export_eset                       = as.logical(default_pm$export_eset)
  env$time_series                       = as.logical(default_pm$time_series) 
  env$quality_control_only              = as.logical(default_pm$quality_control_only)
  env$create_heatmaps_genes_of_interest = as.logical(default_pm$create_heatmaps_genes_of_interest)
  env$use_kegg_for_heatmap              = as.logical(default_pm$use_kegg_for_heatmap)
  env$dif_exp_experiment                = as.logical(default_pm$dif_exp_experiment)
  env$integrate_vcf_files               = as.logical(default_pm$integrate_vcf_files)
  env$multi_probe                       = as.logical(default_pm$multi_probe)
  env$var_filter                        = as.logical(default_pm$var_filter)
  env$integrate_drug_data               = as.logical(default_pm$integrate_drug_data)
  env$use_frma_normalization            = as.logical(default_pm$use_frma_normalization)
  env$heatmap_vis                       = as.logical(default_pm$heatmap_vis)
  env$use_logFc_only                    = as.logical(default_pm$use_logFc_only)
  env$run_generic                       = as.logical(default_pm$run_generic)
  env$annotage_tissue_abbundance        = as.logical(default_pm$annotage_tissue_abbundance)
  env$cohorts_type                      = default_pm$cohorts_type
  env$stat_design                       = default_pm$stat_design
  #absent_genes_file                            = "absent_genes.tab"
  env$res_file                          = default_pm$res_file
  #time_series_res_file                         = "time_series_expression.csv"
  genes_of_interest_file                        = default_pm$genes_of_interest_file
  
  if (! exists("use_kegg_for_heatmap", envir = env)){
    env$use_kegg_for_heatmap = F
  }
  if (! exists("create_heatmaps_genes_of_interest", envir = env)){
    env$create_heatmaps_genes_of_interest = T
  }
  env$heatmap_list_genes_count = 40
  
  # create paths for input and output
  dir.create( paste( output_path, "Output_GeneraPipe", sep = "/"), showWarnings = FALSE )
  cel_files_path                          = sub(x = cel_files_path, "/$", "")
  cpdb_file                               = "CPDB_pathways_genes.tab"
  #tissue_norm_exprs_file                 = "GSE1133-GPL96_series_matrix.txt"
  #tissue_norm_exprs_file_path            = paste( pipeline_loc , paste( "Misc" ,tissue_norm_exprs_file, sep ="/" ) , sep ="/" )
  env$quality_report_path         = paste( output_path, "Output_GeneraPipe/QC_report" , sep = "/" )
  env$results_file_path           = paste( output_path, "Output_GeneraPipe", paste( "Results", project_name, sep ="_"), sep = "/" )
  env$name_res_file               = paste( env$results_file_path, paste( "dif_exp_results", "csv", sep ="."), sep = "/" )
  env$pathway_maps_path           = paste( env$results_file_path  , "pathway_maps_dif_exp", sep ="/" )
  #body_exprs_maps_path                   = paste( pipeline_loc, "Misc/HPM_gene_level_epxression_matrix_Kim_et_al_052914.csv", sep ="/" )
  env$gsea_output_path            = paste( output_path, "Output_GeneraPipe/GSEA_Results", sep = "/" )
  #absent_gene_file_path                  = paste( cel_files_path, absent_genes_file, sep ="/" )
  env$tissue_abbundance_res_file  = paste( cel_files_path, "Tissue_abundance_results.csv", sep ="/" )
  env$kegg_file_path              = paste( cel_files_path , kegg_file, sep ="/")
  env$cpdb_file_path                      = paste(package_path, paste( "extdata", cpdb_file, sep ="/" ), sep ="/")
  #time_series_res_file_path              = paste( cel_files_path, time_series_res_file, sep ="/")
  env$entities_of_interest_path   = paste( env$results_file_path, "Entities_of_interest", sep ="/")
  env$genes_of_interest_file_path = paste( env$entities_of_interest_path, genes_of_interest_file, sep ="/") 
  env$user_folder                 = as.character( system("echo $HOME", intern = T) )
  cohorts_file_path                       = paste( cel_files_path, cohorts_file, sep ="/" )
  
  if (GeneraPipe:::strEndsWith( cohorts_file_path, ".csv")){  
    phenodata = read.table( cohorts_file_path , header = TRUE , sep = "," )
  } else {  
    phenodata = read.table( cohorts_file_path , header = TRUE , sep = "\t" )
  }
  phenodata$ID = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", phenodata$ID )
  
  if (GeneraPipe:::strEndsWith(env$kegg_file_path, ".csv")){  
    env$keggdata = read.table( env$kegg_file_path , header = TRUE , sep = "," )
  } else {  
    env$keggdata = read.table( env$kegg_file_path , header = TRUE , sep = "\t" )
  }
  
  return(phenodata)
}