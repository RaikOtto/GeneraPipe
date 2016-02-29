create_cohorts = function( raw_data, chip_type, set_ctrl, set_case, phenodata){
  
  if (exists("eset")){
  
    if (zipped){
      p_data = paste0( phenodata$ID, ".gz")
    } else {
      p_data = phenodata$ID
    }
  }

  rownames(Biobase::pData(raw_data)) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(Biobase::pData(raw_data))) # the reason is to assure, that no .gz or .CEL ending is present when we add it
  Biobase::sampleNames(raw_data) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", Biobase::sampleNames(raw_data)) 
  #print( paste( "Did not find following sampels in raw data, but in cohorts file: ", paste( !( phenodata$ID[phenodata$ID %in% rownames( pData( raw_data ) ) ] ) , sep =" " ) ) ) 
  phenodata = phenodata[phenodata$ID %in% rownames(Biobase::pData(raw_data)),]
  env$cohorts_vec = as.vector(phenodata[, which(colnames(phenodata) == env$cohorts_type)])

  if (! env$quality_control_only){
  
    env$cohorts_vec[env$cohorts_vec %in% set_ctrl] = "CTRL"
    env$cohorts_vec[env$cohorts_vec %in% set_case] = "CASE"
  
    index_cohorts_vec = which(env$cohorts_vec %in% c("CTRL","CASE") == TRUE)
    env$cohorts_vec = env$cohorts_vec[env$cohorts_vec %in% c("CTRL","CASE")]
    names(env$cohorts_vec) = as.character(phenodata$ID[index_cohorts_vec])
  
    nmbr_samples = sum(env$cohorts_vec %in% c("CTRL","CASE"))
  
    if ( env$stat_design == "intercept" ) {
    
      env$design = model.matrix(~  env$cohorts_vec)  # intercept
    
    }else{
    
      env$design = model.matrix(~ 0 + env$cohorts_vec) # contrast
    }
  
    colnames(env$design)[ colnames( env$design ) == paste0("env$cohorts_vec","CTRL") ] = "CTRL"
    colnames(env$design)[ colnames( env$design ) == paste0("env$cohorts_vec","CASE") ] = "CASE"
  
  }

  if( chip_type == "HumanHT-12.v4" ){
    index_ctrl2 = as.integer( which( phenodata$Group == "NORMAL") )
    index_case2 = as.integer( which( phenodata$Group == "NET" ) )
  }

  env$index_ctrl = as.integer( which( env$design[,colnames(env$design) == "CTRL"] == 1 ) )
  env$index_case = as.integer( which( env$design[,colnames(env$design) == "CASE"] == 1 ) )

  #if (exists("raw_data")){
  if( ! ("Group" %in% colnames(Biobase::pData( raw_data ) ) ) & chip_type != "HumanHT-12.v4" ) {
  
    raw_data_group_vec                             = rep("",dim( Biobase::pData(raw_data) )[1] )
    raw_data_group_vec[env$index_ctrl]     = "CTRL"
    raw_data_group_vec[env$index_case]     = "CASE"
    Biobase::pData(raw_data)                       = cbind(Biobase::pData(raw_data), raw_data_group_vec)
    env$raw_data_group_vec                 = raw_data_group_vec[which(raw_data_group_vec != "")]
    colnames(Biobase::pData(raw_data))[-1]         = "Group"

  } else if( chip_type == "HumanHT-12.v4"){
 
    raw_data_group_vec = rep("",dim( Biobase::pData(raw_data) )[1] )
    raw_data_group_vec[index_ctrl2] = "CTRL"
    raw_data_group_vec[index_case2] = "CASE"
    Biobase::pData(raw_data)$Cohorts = raw_data_group_vec
    env$raw_data_group_vec = raw_data_group_vec[which(raw_data_group_vec != "")]
    eset = raw_data
    eset = eset[ ,c( index_case2, index_ctrl2 )  ]
  }
  
  env$phenodata = phenodata
  return(raw_data)
}

