selected_probe_expression = function(){
  
  d=c("213539_at","205456_at","206804_at","215345_x_at","216298_at","217381_s_at","203434_s_at","203435_s_at","206398_s_at","210356_x_at","213953_at","217418_x_at","228592_at","228599_at","231418_at","206150_at","209543_s_at","205692_s_at","1557905_s_at","1565868_at","204489_s_at","204490_s_at","209835_x_at","210916_s_at","212014_x_at","212063_at","216056_at","217523_at","229221_at","234411_x_at","234418_x_at","1552480_s_at","1569830_at","207238_s_at","212587_s_at","212588_at","209201_x_at","211919_s_at","217028_at")
  
  p_mapping = which( phenodata$ID %in% colnames(eset) )
  p_pheno = phenodata$Group[p_mapping]
  
  f = function(input_vec, grouping ) {
    res = aggregate( input_vec , by = list(grouping), mean );
    print(res[,1])
    return( res[,2] );
  }
  
  grouping_probes_genes = unlist( mget(d, hgu133plus2SYMBOL) )
  grouping_probes_genes[ which( grouping_probes_genes == "TARP")  ]  = "CD3G"
  grouping_probes_genes[ which( grouping_probes_genes == "KRT20")  ] = "CD20"
  grouping_probes_genes[ which( grouping_probes_genes == "MME")  ]   = "CD10"
  grouping_probes_genes[ which( grouping_probes_genes == "MS4A1")  ] = "CD20"
  grouping_probes_genes[ which( grouping_probes_genes == "PTPRC")  ] = "CD45"
  grouping_probes_genes[ which( grouping_probes_genes == "CXCR4")  ] = "CD184"
  
  mapping_prob_row = which( rownames(eset) %in% d)
  res_time_series = exprs(eset)[ mapping_prob_row, ]
  #p_pheno = p_pheno[-excluded_files  ]
  
  groups = c( as.vector( unlist( as.character( aggregate( exprs(eset)[1,] , by = list( p_pheno ), mean )[,1] ) ) ) )
  
  res_mat_probe = cbind( apply( res_time_series,MARGIN=1, FUN = f, p_pheno ) )
  
  res_mat_gene = cbind( apply( t(res_mat_probe) , MARGIN = 2, FUN = f, grouping_probes_genes ) )
  
  colnames(res_mat_gene) = groups
  rownames(res_mat_gene) = c("CD10","CD19","CD20","CD27","CD34","CD38","CD3D","CD3E","CD3G","CD44","CD45","CD184")
  res_mat_gene = round(res_mat_gene,2)
  
  #mapping_order = match(colnames(res_mat_gene), as.vector(mapping_order$Group) )
  res_mat_gene = res_mat_gene[ ,order( mapping_order ) ]
  
  res_mat_gene_final = rbind( colnames(res_mat_gene), res_mat_gene   )
  res_mat_gene_final = cbind(  c( "marker", rownames(res_mat_gene)), res_mat_gene_final)
  
  min_row = apply(res_mat_gene,MARGIN=1,FUN=min)
  rel_res_mat_gene = res_mat_gene - min_row
  max_row = apply(rel_res_mat_gene,MARGIN=1,FUN=max)
  rel_mat = round( rel_res_mat_gene / max_row, 2 )
  
  library(gdata)
  heatmap.2(rel_mat, col=greenred(75), Colv = F)
  
  #res_mat_gene_final = res_mat_gene_final[order( res_mat_gene_final[,1] )  ,]
  write.table( time_series_res_file_path, x =  res_mat_gene_final, row.names=F, col.names=F , quote = F, sep ="\t" )
}
