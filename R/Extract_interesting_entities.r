extract_intr_entities = function( eset, chip_type ){

  dir.create(env$entities_of_interest_path, showWarnings = FALSE)
  #source("Src/cohort_creation.r")

  #if ( ! exists("hgnc_symbols"))
    #source("Src/annotation.r")

  if( chip_type == "HumanHT-12.v4" ){
    expr_data_fix = eset_select
  } else{
    expr_data_fix  = exprs(eset)
  }

  exprs_case = rowMeans(expr_data_fix[, env$index_case] )
  exprs_ctrl = rowMeans(expr_data_fix[, env$index_ctrl] )
  dif_exp    = exprs_case - exprs_ctrl

  expr_data = cbind(
    round(dif_exp, 2),
    round(exprs_case, 2),
    round(exprs_ctrl, 2),
    expr_data_fix
  )

  colnames(expr_data)[1] = "logFC"
  colnames(expr_data)[2] = "Expr_case"
  colnames(expr_data)[3] = "Expr_ctrl"
  expr_data = as.data.frame(expr_data)
  
  cpdb_t = read.table(env$cpdb_file_path, header = TRUE, sep = "\t", fill = TRUE)
  cpdb_ident = stringr::str_replace(cpdb_t$external_id, "path:", "")
  
  interesting_pathways_mapping = match(stringr::str_trim(env$keggdata$Kegg_id), stringr::str_trim(cpdb_ident))
  interesting_pathways_table_kegg_id = cpdb_t[interesting_pathways_mapping, ]
  
  ### genes
  
  if (! is.null(env$keggdata$Gene_id_hgnc)){
    
    selection = as.character(unique(env$keggdata$Gene_id_hgnc))
    selection = selection[selection != ""]
    
    hgnc_symbols = as.character(env$hgnc_genes)
    mapping = which(env$hgnc_genes %in% selection)
    gene_ids = hgnc_symbols[env$hgnc_genes %in% selection]
    
    exprs_case = expr_data_fix[mapping, env$index_case]
    exprs_ctrl = expr_data_fix[mapping, env$index_ctrl]
    dif = rowMeans(exprs_case) - rowMeans(exprs_ctrl)
    
    exprs_case = round(exprs_case, 2)
    exprs_ctrl = round(exprs_ctrl, 2)
    dif = round(dif, 2)
    
    sample_expression = cbind(exprs_ctrl, exprs_case)
    
    if (dim(sample_expression)[2] == 1){
      sample_expression = t(sample_expression)
    }
    
    res_int = cbind(
      as.double(dif),
      as.double(round(rowMeans(exprs_ctrl), 2)),
      as.double(round(rowMeans(exprs_case), 2)),
      sample_expression
    )
    
    colnames(res_int) =  c(
      "logFC",
      "expr_ctrl",
      "expr_case",
      c(
        colnames(eset)[env$index_ctrl],
        colnames(eset)[env$index_case]
      )
    )
    
    res_int = cbind.data.frame(gene_ids, res_int)
    colnames(res_int)[1] = "HGNC_symbol"
    res_int = res_int[order(as.double(res_int[, 2]), decreasing = TRUE ), ]
    
    if (env$integrate_drug_data){
      
      drug_vec = rep( "", length( res_int[1,] )  )
      drug_vec[1] = drug_type
      names( drug_vec ) = colnames( res_int )
      
      if ( ! exists("drug_data") ){
        
        message("Did not find drug data, but drug_response parameter was set")
        exit()
      }
      
      mapping = match( names( drug_vec ), names(drug_data), nomatch = 0  )
      
      if ( ! (T %in% ( mapping != 0 ) ) ){
        
        message("Problem mapping samplenames of genes of interest to drug data from cohorts file")
        exit()
      }
      
      drug_vec[ which( mapping != 0  ) ] = as.double( drug_data[ mapping  ] )
      res_ctrl  = as.double( drug_data[ index_ctrl ] )
      res_case = as.double( drug_data[ index_case ] )
      dif = round( mean( res_case ) - mean(as.double( res_ctrl  )), 1 )
      
      # drug response
      
      drug_vec[ match( "expr_case", names(drug_vec) ) ] = mean( res_case )
      drug_vec[ match( "expr_ctrl", names(drug_vec) ) ] = mean( res_ctrl )
      drug_vec[ match( "logFC", names(drug_vec) ) ]     = dif
      
      # correlations
      index = seq( 5, dim(res_int)[2]  )
      cor_pearson  = c()
      cor_spearman = c()
      
      for (n in seq( 1,dim(res_int)[1]  ) ){
        
        cor_pearson  = c( cor_pearson, cor( as.double( res_int[n,  index ] ), as.double( drug_vec[index] ), method = "pearson", use = "pairwise.complete.obs"))
        cor_spearman = c( cor_spearman, cor( as.double( res_int[n,  index ] ), as.double( drug_vec[index] ), method = "spearman", use = "pairwise.complete.obs"))
      }
      
      cor_pearson = round(cor_pearson, 2)
      cor_spearman = round(cor_spearman, 2)
      
      # cohorts
      
      case_ctrl = phenodata [ match( names(cohorts_vec), phenodata$ID  ) , which( colnames(phenodata) == cohorts_type  )  ]
      case_ctrl = c( "Cohort", "" , "" , "", as.character(case_ctrl)  )
      
      res_int = rbind( res_int, drug_vec  )
      res_int = rbind( res_int, case_ctrl  )
      
      res_int = cbind( res_int[ ,1:4], cor_pearson, cor_spearman, res_int[ , 5:dim(res_int)[2]  ] )
    }
    
    xlsx::write.xlsx(res_int, stringr::str_replace(stringr::str_replace(env$genes_of_interest_file_path, "~", env$user_folder), ".csv", ".xls"), row.names = FALSE)
  }
  
  ### pathways
  
  if (! is.null(env$keggdata$Kegg_id)){
    
    mapping = match(as.character(env$keggdata$Kegg_id), cpdb_ident, nomatch = 0 )
    message(paste0(c("Not machted Kegg Pathways:", as.character(env$keggdata$Kegg_id[mapping == 0])), collapse = ", "))
    
    #mapping = match( cpdb_ident, keggdata$Kegg_id, nomatch = 0 )
    #mapping = mapping[mapping!=0]
    genes_of_interest = cpdb_t$hgnc_symbol_ids[mapping]
    
    for (i in mapping){
      
      pathway_id = as.character(cpdb_ident[0])
      
      if (length(pathway_id) >= 1){
        
        pathway_name  = as.character(pathway_id)
        pathway_name  = stringr::str_replace(pathway_name, "/", "and")
        genes         = cpdb_t$hgnc_symbol_ids[i]
        gene_list     = unlist(stringr::str_split(genes, ","))
        mapping_gene  = which(as.character(hgnc_symbols) %in% as.character(gene_list))  
        exprs_genes   = expr_data_fix[mapping_gene, ]
        exprs_case    = round(rowMeans(exprs_genes[, env$index_case]), 2)
        exprs_ctrl    = round(rowMeans(exprs_genes[, env$index_ctrl]), 2)
        dif_exp       = round(exprs_case - exprs_ctrl, 2)
        
        exprs_genes = cbind.data.frame(as.double(dif_exp), as.double(exprs_ctrl), as.double(exprs_case), hgnc_symbols[mapping_gene] , round(expr_data_fix[mapping_gene, ], 2))
        colnames(exprs_genes) =  c("logFC", "expr_ctrl", "expr_case", "hgnc_symbol", colnames(expr_data_fix))
        exprs_genes = exprs_genes[order(as.double(exprs_genes[, 1]), decreasing = TRUE), ]
        
        if (env$integrate_drug_data ){
          drug_vec = rep( "", length( exprs_genes[1,] )  )
          drug_vec[1] = drug_type
          names( drug_vec ) = colnames( exprs_genes )
          
          if ( ! exists("drug_data") ){
            
            message("Did not find drug data, but drug_response parameter was set")
            exit()
          }
          
          mapping = match( names( drug_vec ), names(drug_data), nomatch = 0  )
          
          if ( ! (T %in% ( mapping != 0 ) ) ){
            
            message("Problem mapping samplenames of genes of interest to drug data from cohorts file")
            exit()
          }
          
          drug_vec[ which( mapping != 0  ) ] = as.double( drug_data[ mapping  ] )
          res_ctrl  = as.double( drug_data[ index_ctrl ] )
          res_case = as.double( drug_data[ index_case ] )
          dif = round( mean( res_case ) - mean(as.double( res_ctrl  )), 1 )
          
          # drugs
          
          drug_vec[ match( "expr_case", names(drug_vec) ) ] = mean( res_case )
          drug_vec[ match( "expr_ctrl", names(drug_vec) ) ] = mean( res_ctrl )
          drug_vec[ match( "logFC", names(drug_vec) ) ]     = dif
          drug_vec[ match( "hgnc_symbol", names(drug_vec) ) ] = drug_type
          
          # corrs
          
          index = seq( 5, dim(exprs_genes)[2]  )
          cor_pearson  = c()
          cor_spearman = c()
          
          for (n in seq( 1,dim(exprs_genes)[1]  ) ){
            
            cor_pearson  = c( cor_pearson, cor( as.double( exprs_genes[n,  index ] ), as.double( drug_vec[index] ), method = "pearson", use = "pairwise.complete.obs"))
            cor_spearman = c( cor_spearman, cor( as.double( exprs_genes[n,  index ] ), as.double( drug_vec[index] ), method = "spearman", use = "pairwise.complete.obs"))
          }
          
          cor_pearson = round(cor_pearson, 2)
          cor_spearman = round(cor_spearman, 2)
          
          # cohorts
          
          case_ctrl = phenodata [, which( colnames(phenodata) == cohorts_type  )  ]
          case_ctrl = c( "Cohort", "" , "" , "", as.character(case_ctrl)  )
          
          exprs_genes = rbind( exprs_genes, drug_vec  )
          exprs_genes = rbind( exprs_genes, case_ctrl  )
          
          exprs_genes = cbind( exprs_genes[ ,1:4], cor_pearson, cor_spearman, exprs_genes[ , 5:dim(exprs_genes)[2]  ] )
        }
        
        file_name = stringr::str_replace(env$genes_of_interest_file_path, "genes_of_interest", paste(pathway_id, pathway_name, sep = "_"))
        message(c(i, file_name))
        
        xlsx::write.xlsx(exprs_genes, stringr::str_replace(stringr::str_replace(file_name, "~", env$user_folder), ".csv", ".xls"), row.names = FALSE)
        
      }
    }
  }
  
  #if ( create_heatmaps_genes_of_interest  ){
    
  #  source("Src/create_heatmaps.r")
  #}
}
  