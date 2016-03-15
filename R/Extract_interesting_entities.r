extract_intr_entities = function(results, paths, pheno_kegg, chip_type, annotation){

  eset = results@eset
  index_case = results@index_case
  index_ctrl = results@index_ctrl
  keggdata = pheno_kegg$kegg
  entities_of_interest_path = paths@entities_of_interest_path
  cpdb_file_path = paths@cpdb_file_path
  genes_of_interest_file_path = paths@genes_of_interest_file_path
  user_folder = paths@user_folder
  hgnc_genes = annotation@hgnc_genes
  
  dir.create(entities_of_interest_path, showWarnings = FALSE)

  exprData  = exprs(eset)
  exprs_case = rowMeans(exprData[, index_case] )
  exprs_ctrl = rowMeans(exprData[, index_ctrl] )
  dif_exp    = exprs_case - exprs_ctrl

  expr_data = cbind(
    round(dif_exp, 2),
    round(exprs_case, 2),
    round(exprs_ctrl, 2),
    exprData
  )

  colnames(expr_data)[1] = "logFC"
  colnames(expr_data)[2] = "Expr_case"
  colnames(expr_data)[3] = "Expr_ctrl"
  expr_data = as.data.frame(expr_data)
  
  cpdb_t = read.table(cpdb_file_path, header = TRUE, sep = "\t", fill = TRUE)
  cpdb_ident = stringr::str_replace(cpdb_t$external_id, "path:", "")
  
  interesting_pathways_mapping = match(stringr::str_trim(keggdata$Kegg_id), stringr::str_trim(cpdb_ident))
  interesting_pathways_table_kegg_id = cpdb_t[interesting_pathways_mapping, ]
  
  ### genes
  
  if (! is.null(keggdata$Gene_id_hgnc)){
    
    selection = as.character(unique(keggdata$Gene_id_hgnc))
    selection = selection[selection != ""]
    
    hgnc_symbols = as.character(hgnc_genes)
    mapping = which(hgnc_genes %in% selection)
    gene_ids = hgnc_symbols[hgnc_genes %in% selection]
    
    exprs_case = exprData[mapping, index_case]
    exprs_ctrl = exprData[mapping, index_ctrl]
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
        colnames(eset)[index_ctrl],
        colnames(eset)[index_case]
      )
    )
    
    res_int = cbind.data.frame(gene_ids, res_int)
    colnames(res_int)[1] = "HGNC_symbol"
    res_int = res_int[order(as.double(res_int[, 2]), decreasing = TRUE ), ]
    
    xlsx::write.xlsx(res_int, stringr::str_replace(stringr::str_replace(genes_of_interest_file_path, "~", user_folder), ".csv", ".xls"), row.names = FALSE)
  }
  
  ### pathways
  
  if (! is.null(keggdata$Kegg_id)){
    
    mapping = match(as.character(keggdata$Kegg_id), cpdb_ident, nomatch = 0 )
    message(paste0(c("Not machted Kegg Pathways:", as.character(keggdata$Kegg_id[mapping == 0])), collapse = ", "))
    
    #mapping = match( cpdb_ident, keggdata$Kegg_id, nomatch = 0 )
    #mapping = mapping[mapping!=0]
    genes_of_interest = cpdb_t$hgnc_symbol_ids[mapping]
    
    for (i in mapping){
      
      pathway_id = as.character(cpdb_ident[i])
      pathway_name = as.character(cpdb_t$pathway[i])
      
      if (length(pathway_id) >= 1){
        
        #pathway_name  = as.character(pathway_id)
        pathway_name  = stringr::str_replace(pathway_name, "/", "and")
        genes         = cpdb_t$hgnc_symbol_ids[i]
        gene_list     = unlist(stringr::str_split(genes, ","))
        mapping_gene  = which(as.character(hgnc_symbols) %in% as.character(gene_list))  
        exprs_genes   = exprData[mapping_gene, ]
        exprs_case    = round(rowMeans(exprs_genes[, index_case]), 2)
        exprs_ctrl    = round(rowMeans(exprs_genes[, index_ctrl]), 2)
        dif_exp       = round(exprs_case - exprs_ctrl, 2)
        
        exprs_genes = cbind.data.frame(as.double(dif_exp), as.double(exprs_ctrl), as.double(exprs_case), hgnc_symbols[mapping_gene] , round(exprData[mapping_gene, ], 2))
        colnames(exprs_genes) =  c("logFC", "expr_ctrl", "expr_case", "hgnc_symbol", colnames(exprData))
        exprs_genes = exprs_genes[order(as.double(exprs_genes[, 1]), decreasing = TRUE), ]
        
        file_name = stringr::str_replace(genes_of_interest_file_path, "genes_of_interest", paste(pathway_id, pathway_name, sep = "_"))
        #message(c(i, file_name))
        
        xlsx::write.xlsx(exprs_genes, stringr::str_replace(stringr::str_replace(file_name, "~", user_folder), ".csv", ".xls"), row.names = FALSE)
      }
    }
  }
  
}
  