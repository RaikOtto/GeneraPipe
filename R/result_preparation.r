result_preparation = function(paths, results, annotation, chip_type, lfc_exp){
  
  topall        = results@topall
  volc_all      = results@volc_all
  design        = results@design
  eset          = results@eset
  index_case    = results@index_case
  index_ctrl    = results@index_ctrl
  index_topall  = results@index_topall
  hgnc_genes    = annotation@hgnc_genes
  hgnc_names    = annotation@hgnc_names
  ensembl_genes = annotation@ensembl_genes
  
  #### generic output
  #probe_ids = rownames(topall)
  
  #index_probes = match(rownames(topall), rownames(eset), nomatch = 0)
  #exprs_case = rowMeans(exprs(eset)[index_probes, index_case])
  #exprs_ctrl = rowMeans(exprs(eset)[index_probes, index_ctrl])
  
  #hgnc_genes    = hgnc_genes[index_topall]
  #hgnc_names    = hgnc_names[index_topall]
  #ensembl_genes = ensembl_genes[index_topall]
  
  #topall_res = data.frame(
  #  "Probe_ids"           = probe_ids,
  #  "expr_ctrl"           = round(exprs_ctrl, 2),
  #  "expr_case"           = round(exprs_case, 2),
  #  "logFC"               = topall$logFC,
  #  "P_Value"             = topall$P.Val,
  #  "HGNC_symb"           = hgnc_genes,
  #  "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";")
  #)
  #topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE), ]
  
  probe_ids    = rownames(volc_all)
  index_probes = match(rownames(volc_all), rownames(eset), nomatch = 0)
  exprs_case   = rowMeans(exprs(eset)[index_probes, index_case])
  exprs_ctrl   = rowMeans(exprs(eset)[index_probes, index_ctrl])
  
  hgnc_genes = hgnc_genes[index_probes]
  hgnc_names = hgnc_names[index_probes]
  ensembl_genes = unlist(ensembl_genes[index_probes])
  
  topall_res = data.frame(
    "Probe_ids"           = probe_ids,
    "expr_ctrl"           = round(exprs_ctrl, 2),
    "expr_case"           = round(exprs_case, 2),
    "logFC"               = round(volc_all$logFC,2),
    "P_Value"             = round(volc_all$adj.P.Val, 4),
    "HGNC_symb"           = hgnc_genes,
    "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";")
  )
  
  topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE), ]
  
  if(FALSE){
  ## unify if possible output for all different chip types
  if (chip_type == "hgu133plus2"){
    probe_ids    = row.names(topall)
    hgnc_genes   = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2SYMBOL))
    hgnc_names   = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2GENENAME))
    entrez_genes = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2ENTREZID))
    pathway      = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2PATH))
    index_probes = match(rownames(topall), rownames(eset), nomatch = 0)
    exprs_case   = rowMeans(exprs(eset)[index_probes, index_case])
    exprs_ctrl   = rowMeans(exprs(eset)[index_probes, index_ctrl])

    topall_res = data.frame(
      "Probe_ids"           = probe_ids,
      "expr_ctrl"           = round(exprs_ctrl, 2),
      "expr_case"           = round(exprs_case, 2),
      "logFC"               = topall$logFC,
      "P_Value"             = topall$P.Val,
      "HGNC_symb"           = hgnc_genes,
      "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";"),
      "pathway"             = pathway
    )
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE), ]

  } else if (chip_type == "hgu133a") {
    hgnc_genes   = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aSYMBOL))
    hgnc_names   = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aGENENAME))
    entrez_genes = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aENTREZID))
    pathway      = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aPATH))
    
    index_probes = match(rownames(topall), rownames(eset), nomatch = 0)
    exprs_case   = rowMeans(exprs(eset)[index_probes, index_case])
    exprs_ctrl   = rowMeans(exprs(eset)[index_probes, index_ctrl])

    topall_res = data.frame(
      "logFC"               = topall$logFC,
      "expr_ctrl"           = round(exprs_ctrl, 2),
      "expr_case"           = round(exprs_case, 2),
      "P_Value"             = topall$P.Val,
      "HGNC_symb"           = hgnc_genes,
      "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";"),
      "pathway"             = pathway
    )
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]

  } else if (chip_type %in% c("pd.hugene.2.0.st", "pd.huex.1.0.st.v2")){
    probe_ids    = rownames(topall)
    index_probes = match(rownames(topall), rownames(eset), nomatch = 0)
    exprs_case   = rowMeans(exprs(eset)[index_probes, index_case])
    exprs_ctrl   = rowMeans(exprs(eset)[index_probes, index_ctrl])
    hgnc_symbols = stringr::str_trim(unlist(lapply(topall$geneassignment, FUN = GeneraPipe:::split_fun, 2)))
    hgnc_names   = stringr::str_trim(unlist(lapply(topall$geneassignment, FUN = GeneraPipe:::split_fun, 3)))

    topall_res = data.frame(
      "Probe_ids"           = probe_ids,
      "logFC"               = round(topall$logFC, 2),
      "expr_ctrl"           = round(exprs_ctrl, 2),
      "expr_case"           = round(exprs_case, 2),
      "P_Value"             = topall$P.Value,
      "HGNC_symb"           = hgnc_symbols,
      "HGNC_names"          = hgnc_names,
      "Gene_assignment"     = topall$geneassignment
    )
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]
  } 
  }

  dir.create(paths@results_file_path, showWarnings = FALSE)
  xlsx::write.xlsx(topall_res, stringr::str_replace(stringr::str_replace(paths@name_res_file, "~", paths@user_folder),".csv", ".xls"), row.names = FALSE)
  
  results@topall_res = topall_res
  return(results)
}