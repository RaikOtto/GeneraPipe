create_cohorts = function( results, chip_type, set_ctrl, set_case, pheno_kegg, stat_design){
  
  raw_data = results@raw_data
  phenodata = pheno_kegg$pheno
  
  rownames(Biobase::pData(raw_data)) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", rownames(Biobase::pData(raw_data))) # the reason is to assure, that no .gz or .CEL ending is present when we add it
  Biobase::sampleNames(raw_data) = base::gsub(c(".gz|.CEL|.cel|.GZ"), "", Biobase::sampleNames(raw_data)) 
  phenodata = phenodata[phenodata$ID %in% rownames(Biobase::pData(raw_data)),]
  cohorts_vec = as.vector(phenodata[, which(colnames(phenodata) == "Group")])
  
  cohorts_vec[cohorts_vec %in% set_ctrl] = "CTRL"
  cohorts_vec[cohorts_vec %in% set_case] = "CASE"
  
  index_cohorts_vec = which(cohorts_vec %in% c("CTRL", "CASE") == TRUE)
  cohorts_vec = cohorts_vec[cohorts_vec %in% c("CTRL", "CASE")]
  names(cohorts_vec) = as.character(phenodata$ID[index_cohorts_vec])
  
  nmbr_samples = sum(cohorts_vec %in% c("CTRL", "CASE"))
  
  if (stat_design == "intercept") {
    
    design = model.matrix(~ env$cohorts_vec)  # intercept
    
  }else{
    
    design = model.matrix(~ 0 + cohorts_vec) # contrast
  }
  
  colnames(design)[colnames(design) == paste0("cohorts_vec","CTRL") ] = "CTRL"
  colnames(design)[colnames(design) == paste0("cohorts_vec","CASE") ] = "CASE"

  index_ctrl = as.integer(which(design[, colnames(design) == "CTRL"] == 1))
  index_case = as.integer(which(design[, colnames(design) == "CASE"] == 1))

  
  if(! ("Cohort" %in% colnames(Biobase::pData(raw_data))) & chip_type != "HumanHT-12.v4") {
    raw_data_group_vec                     = rep("", dim(Biobase::pData(raw_data))[1])
    raw_data_group_vec[index_ctrl]     = "CTRL"
    raw_data_group_vec[index_case]     = "CASE"
    Biobase::pData(raw_data)               = cbind(Biobase::pData(raw_data), raw_data_group_vec)
    raw_data_group_vec                     = raw_data_group_vec[which(raw_data_group_vec != "")]
    colnames(Biobase::pData(raw_data))[-1] = "Cohort"
  }

  results@raw_data    = raw_data
  results@cohorts_vec = cohorts_vec
  results@design      = design
  results@index_case  = index_case
  results@index_ctrl  = index_ctrl
  results
}