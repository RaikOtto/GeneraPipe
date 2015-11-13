
run_analysis = function(){
  
  if ( run_generic ){
    
    message( "Running step 1: Generic parameter init" )
    source( "Src/set_generic_initial_parameters.r" )
  }
  
  ### Part Parse Cel Files
  
  if ( parse_files ){
    
    message( "Running step 2: Parsing raw CEL files")
    source("Src/parse_cel_files.r")
  }
  
  ### Part cohort creation
  
  if ( create_cohorts ){
    
    message( "Running step 3: Cohort construction" )
    source("Src/cohort_creation.r")
  }
  
  ### Part Normalization
  
  if ( normalize ){
    if ( use_frma_normalization ){
      source("Src/frma_normalization.r")
    }else{
      
      message( "Running step 4: Normalizing" )
      source("Src/normalization.r")
    }
  }
  
  # Quality Control
  
  if ( qc_control ){
    
    message( "Running step 5: Quality Control")
    source("Src/quality_control.r")
  }
  
  ### Part annotation
  
  if ( annotate ){
    
    message( "Running step 6: Data annotation" )
    source("Src/annotation.r")
  }
  
  ### Part Excluding non-present Genes
  
  if (absent_analysis){
    
    message( "Running step 7: Absent gene analysis" )
    source("Src/absent_genes_analysis.r")
  }
  
  ### Part Differential Expression
  
  if ( dif_exp_ana ){
    
    message( "Running step 8: Differential Expression" )
    if ( use_gsea ){
      source( "Src/gsea.r")
    } else{
      source("Src/differential_expression.r")
    }
  }
  
  ### Part Analysis
  
  if ( heatmap_vis  )
    source( "Src/heatmap_probes_of_interest.r" )
  
  ### Part Visualization
  
  if ( export_results ){
    
    message( "Running step 9: Exporting data" )
    source( "Src/result_preparation.r" )
  }
  
  if ( create_pathways ){
    
    if (! time_series  ){
      
      message( "Running step 10: Creating Pathway maps" )
      source("Src/create_pathway_maps.r")
    }
  }
  
  if (extract_interest){
    
    message( "Running step 11: Extraction of interesting information" )
    source("Src/Extract_interesting_entities.r")
  }
  
  if (annotate_tissue_abbundance){
    
    message( "Comparison to normal tissue expression" )
    source("Src/annotate_tissue_abbundance.r")
  }
}