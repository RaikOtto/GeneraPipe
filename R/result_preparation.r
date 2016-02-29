result_preparation = function( eset, topall, chip_type, lfc_exp ){
  
  if (chip_type == "hgu133plus2"){
    probe_ids    = row.names(topall)
    hgnc_genes   = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2SYMBOL))
    hgnc_names   = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2GENENAME))
    entrez_genes = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2ENTREZID))
    pathway      = as.character(BiocGenerics::mget(rownames(topall), hgu133plus2.db::hgu133plus2PATH))

    topall_res = data.frame(
      "logFC"               = topall$logFC,
      "P_Value"             = topall$P.Val,
      "HGNC_symb"           = hgnc_genes,
      "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";"),
      "entrez"              = entrez_genes,
      "pathway"             = pathway
    )

    row.names(topall_res) = probe_ids
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]

  } else if (chip_type == "hgu133a") {
    hgnc_genes   = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aSYMBOL))
    hgnc_names   = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aGENENAME))
    entrez_genes = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aENTREZID))
    pathway      = as.character(BiocGenerics::mget(rownames(topall), hgu133a.db::hgu133aPATH))

    topall_res = data.frame(
      "logFC"               = topall$logFC,
      "P_Value"             = topall$P.Val,
      "HGNC_symb"           = hgnc_genes,
      "HGNC_name"           = stringr::str_replace_all(hgnc_names,",",";"),
      "entrez"              = entrez_genes,
      "pathway"             = pathway
    )

    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]

  } else if (chip_type %in% c("pd.hugene.2.0.st", "pd.huex.1.0.st.v2")){
   # if ( ! exists("index_case"))
      #source("Src/annotation.r")

    probe_ids    = rownames(topall)
    index_probes = match(rownames(topall), rownames(eset), nomatch = 0)
    exprs_case   = rowMeans(exprs(eset)[index_probes, env$index_case])
    exprs_ctrl   = rowMeans(exprs(eset)[index_probes, env$index_ctrl])
    hgnc_ids     = stringr::str_trim(unlist(lapply(topall$geneassignment, FUN = GeneraPipe::split_fun, 2)))
    hgnc_names   = stringr::str_trim(unlist(lapply(topall$geneassignment, FUN = GeneraPipe::split_fun, 3)))

    topall_res = data.frame(
      "logFC"               = round(topall$logFC, 2),
      "expr_ctrl"           = round(exprs_ctrl, 2),
      "expr_case"           = round(exprs_case, 2),
      "P_Value"             = topall$P.Value,
      "HGNC_symb"           = hgnc_ids,
      "HGNC_names"          = hgnc_names,
      "Gene_assignment"     = topall$geneassignment,
      "Probe_ids"           = probe_ids
    )

    row.names(topall_res) = probe_ids
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]

  } else if ( chip_type == "HumanHT-12.v4" ){
    gpl = annotation(eDatSet)
    platf = getGEO(gpl, AnnotGPL = TRUE)
    ncbifd = data.frame(attr(dataTable(platf), "table"))

    topall_res = topall[setdiff(colnames(topall), setdiff(fvarLabels(eDatSet), "ID"))]
    topall_res = merge(topall_res, ncbifd, by = "ID")
    topall_res = topall_res[order(topall_res$logFC, decreasing = TRUE),]
  }

  dir.create(env$results_file_path, showWarnings = FALSE)
  xlsx::write.xlsx(topall_res, stringr::str_replace(stringr::str_replace(env$name_res_file, "~", env$user_folder),".csv", ".xls"), row.names = FALSE)
  
  message(c("Amount genes higher in Case cohort: ", sum(topall_res$logFC >= lfc_exp)))
  message(c("Amount genes lower in Case cohort: " , sum(topall_res$logFC <= lfc_exp)))
  
  return(topall_res)
}