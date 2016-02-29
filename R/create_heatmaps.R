create_heatmaps = function( eset, topall_res, project_name, set_ctrl, set_case ){
  
  if (env$create_heatmaps_genes_of_interest){
    
    if (env$use_kegg_for_heatmap){
      
      gene_source_kegg = as.character(env$keggdata$Heatmap[env$keggdata$Heatmap != ""])
      mapping   = which(as.character(env$hgnc_symbols) %in% as.character(gene_source_kegg))
      probe_ids = rownames(exprs(eset))[mapping]
      eset_selection = exprs(eset)[match(probe_ids, rownames(exprs(eset))), ]
      rownames(eset_selection) = env$hgnc_symbols[mapping]
      
    } else {
      
      cutoff = sort(abs(topall_res$logFC), decreasing = TRUE)[env$heatmap_list_genes_count]
      
      probe_ids = rownames(topall_res)[abs(topall_res$logFC) >= cutoff]
      hgnc_ids  = topall_res$HGNC_symb[abs(topall_res$logFC) >= cutoff]
      eset_selection = exprs(eset)[match(probe_ids, rownames(exprs(eset))), ]
      BiocGenerics::rownames(eset_selection) = hgnc_ids
    }
    
    # mapping to cohorts vectors
    
    map_case_ctrl = match(names(env$cohorts_vec), colnames(eset_selection))
    index_ctrl    = which(colnames(eset_selection) %in% names(env$cohorts_vec[env$cohorts_vec == "CTRL"]))
    index_case    = which( colnames(eset_selection) %in% names(env$cohorts_vec[env$cohorts_vec == "CASE"]))
    
    eset_selection = eset_selection[rownames(eset_selection) != "NA", ]
    
    dif = as.double(rowMeans(eset_selection[, index_case]) - rowMeans(eset_selection[, index_ctrl]))
    
    eset_selection = eset_selection[ 
      order(dif, decreasing = TRUE),
      order(order(c(index_ctrl, index_case)))
      ]
    eset_selection_dif = eset_selection - rowMeans(eset_selection)
    
    ## filter for uniques
    
    unique_mapping = match(unique(rownames(eset_selection)), rownames(eset_selection))
    eset_selection = eset_selection[unique_mapping, ]
    eset_selection_dif = eset_selection_dif[unique_mapping, ]
    dif = dif[unique_mapping]
    
    ## trim
    
    for (i in 1:dim(eset_selection_dif)[1]){
      for (j in 1:dim(eset_selection_dif)[2]){
        if (eset_selection_dif[i,j] > 5)
          eset_selection_dif[i,j] = 5
        else if (eset_selection_dif[i,j] < -5)
          eset_selection_dif[i,j] = -5
      }
    }
    
    pdf_name = paste(
      output_path,
      paste0(
        paste0(
          "Output_GeneraPipe/",
          paste0("Results_", project_name)
        ),
        "/heatmap.pdf"
      ),
      sep = "/"
    )
    
    # heatmap.3
    
    subtype_top_bar = as.matrix(c(
      grDevices::colorRampPalette(colors = c("blue"))(length(map_case_ctrl))
    ))
    subtype_top_bar[which(colnames(eset_selection) %in% names(env$cohorts_vec[env$cohorts_vec == "CASE" ]) ), 1  ] = "yellow"
    colnames(subtype_top_bar) = c("Cohort")
    
    logFC_side_bar = t(
      c(
        grDevices::colorRampPalette(colors = c("red"))( length( dif[ dif > 0 ]) ),
        grDevices::colorRampPalette(colors = c("yellow"))( length( dif[ dif == 0 ]) ),
        grDevices::colorRampPalette(colors = c("green"))( length( dif[ dif < 0 ]) )
      )
    )
    rownames(logFC_side_bar) = c("FoldChange")
    library("devtools")
    devtools::source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
    m = grDevices::colorRampPalette(colors = c("green","black","red"))(75)
    
    ### dif darstellung
    
    pdf_name = paste(
      output_path,
      paste0(
        paste0(
          "Output_GeneraPipe/",
          paste0( "Results_", project_name)
        ),
        "/heatmap_difference_overall_expression.pdf"
      ),
      sep = "/"
    )
    
    # heatmap.3
    
    pdf(pdf_name)
    
    heatmap.3(
      eset_selection_dif,
      main = paste0("Sample exprs and logFC ", project_name),
      RowSideColors = logFC_side_bar,
      col = m,
      trace = NULL,
      Colv = FALSE,
      Rowv = FALSE,
      dendrogram = "none" ,
      #margins = c(9,7),
      cexCol = .75,
      cexRow = .75,
      ColSideColors = subtype_top_bar
      
    )
    legend("topright",
           legend=c(
             #"CTRL",
             paste0(c("Control cohort (", set_ctrl,")"), collapse = ""),
             #"CASE",
             paste0(c( "Case cohort (", set_case ,")"), collapse = ""),
             "",
             #"Stronger exp CASE",
             paste0(c("Stronger exp case (", set_case ,")"), collapse = ""),
             "No differential exp",
             #"Stronger exp CTRL"
             paste0( c( "Stronger exp control (", set_ctrl ,")"), collapse = "")
           ),
           fill = c("blue", "yellow", "white", "red", "yellow", "green"), 
           border = FALSE, bty = "n",
           y.intersp = 0.7,
           cex = 0.7
    )
    dev.off()
  }
}