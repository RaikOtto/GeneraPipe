create_pathway_maps = function( topall_res, eset, volc_all, chip_type, package_path ){

  if (! dir.exists(env$pathway_maps_path)){
    dir.create(env$pathway_maps_path)
  } 
  #if (! exists("eset")){
    #source("Src/normalization.r")}

  if (chip_type %in% c("hgu133plus2")){
    hgnc_symbols = as.character(unlist(BiocGenerics::mget(rownames(volc_all), hgu133plus2.db::hgu133plus2SYMBOL)))
  
  } else if ( chip_type %in% c("hgu133a")){
    hgnc_symbols = as.character(unlist(BiocGenerics::mget(rownames(volc_all), hgu133a.db::hgu133aSYMBOL)))
  
  } else if (chip_type %in% c("HumanHT-12.v4")){
    hgnc_symbols = ncbifd$Gene.symbol
    hgnc_symbols = hgnc_symbols[ 93:length(hgnc_symbols)]
  
  } else if ( chip_type %in% c("pd.hugene.2.0.st", "pd.huex.1.0.st.v2")){
    if ((! ("HGNC_symb" %in% colnames(topall_res))) | (length(topall_res$HGNC_symb) == 0)){
      print("Pathway map creation error: topall_res$HGNC_symb does no exist")
      quit()
    }
    hgnc_symbols = stringr::str_trim(unlist(lapply(volc_all$geneassignment, FUN = GeneraPipe:::split_fun, 2)))
    hgnc_names   = stringr::str_trim(unlist(lapply(volc_all$geneassignment, FUN = GeneraPipe:::split_fun, 3)))

  } else{ 
    # note: check if this works for affy 
    hgnc_symbols = stringr::str_trim(unlist(lapply(fData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 2)))
    hgnc_names   = stringr::str_trim(unlist(lapply(fData(eset)$geneassignment, FUN = GeneraPipe:::split_fun, 3)))  
  }

  ensembl     = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  entrez_ids  = biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), values = unique(hgnc_symbols), filters = "hgnc_symbol", mart = ensembl)
  hgnc_map    = match(hgnc_symbols, entrez_ids$hgnc_symbol, nomatch = 0)
  entrez      = rep("", length(hgnc_symbols))
  entrez[which(hgnc_map != 0)] = entrez_ids$entrezgene[hgnc_map]
  entrez[is.na(entrez)] = ""

  if (chip_type == "HumanHT-12.v4"){
    probe_ids   = ncbifd$ID
    eset_select = exprs(eset)[match(probe_ids, rownames(exprs(eset))), ]
    eset_select = eset_select[which( !is.na( rownames(eset_select) ) ), ]
    exprs_case  = round(rowMeans(eset_select[, env$index_case]), 2)
    exprs_ctrl  = round(rowMeans(eset_select[, env$index_ctrl]), 2)
    dif_exp     = round(exprs_case - exprs_ctrl, 2)
  
  } else{
    exprs_case = round(rowMeans(exprs(eset)[, env$index_case]), 2 )
    exprs_ctrl = round(rowMeans(exprs(eset)[, env$index_ctrl]), 2 )
    dif_exp    = round(exprs_case - exprs_ctrl, 2)
  }

  res_all_path = data.frame(
    "logFC"  = dif_exp,
    "HGNC"   = hgnc_symbols,
    "entrez" = entrez
  )

  setwd(env$pathway_maps_path)
  Kegg_id = as.character(env$keggdata$Kegg_id)

  ent = res_all_path$entrez
  ent = ent[ent != ""]
  exp_dat = res_all_path$logFC[res_all_path$entrez != ""]
  names(exp_dat) = res_all_path$entrez[res_all_path$entrez != ""]

  kegg = BiocGenerics::mget(as.character(ent), KEGG.db::KEGGEXTID2PATHID, ifnotfound = list(NA))

  library("pathview")
  for (id in Kegg_id){
    pathview::pathview(gene.data = exp_dat, pathway.id = id, kegg.native = TRUE, limit = max(abs(exp_dat)))
  }

  system("rm *.xml")
  system("ls | grep -v 'pathview.png$' | xargs rm")

  setwd(package_path)
}