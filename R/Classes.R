setClass("project",
  representation(
    ID           = "character",
    project_name = "character",
    chip_type    = "character",
    cohorts_file = "character",
    set_ctrl     = "character",
    set_case     = "character",
    kegg_file    = "character",
    zipped       = "logical",
    p_val        = "numeric",
    lfc_exp      = "numeric",
    stat_design  = "character"
    ),

  prototype(
    ID           = "GSE14407",
    project_name = "Ovarian_GSE14407_test",
    chip_type    = "hgu133plus2",
    cohorts_file = "cohorts.csv",
    set_ctrl     = "Normal",
    set_case     = "Tumor",
    kegg_file    = "pathways.csv",
    zipped       = TRUE, 
    p_val        = .001, 
    lfc_exp      = 1,
    stat_design  = "contrast"
    )
)

setClass("paths",
  representation(
    quality_report_path = "character",
    results_file_path   = "character",
    name_res_file       = "character",
    pathway_maps_path   = "character",
    gsea_output_path    = "character",
    kegg_file_path      = "character",
    cpdb_file_path      = "character",
    entities_of_interest_path   = "character",
    genes_of_interest_file_path = "character",
    user_folder = "character",
    cohorts_file_path = "character",
    package_path = "character",
    gene_set_db_path = "character")
)


setClass("results",
  representation(
    raw_data = "ExonFeatureSet",
    eset = "ExpressionSet",
    topall = "data.frame",
    volc_all = "data.frame",
    topall_res = "data.frame",
    cohorts_vec = "character",
    design = "matrix",
    index_case = "integer",
    index_ctrl = "integer")
)

setClass("annotation",
  representation(
    hgnc_genes        = "character",
    hgnc_symbols  = "character",
    hgnc_names = "character",
    ensembl_genes = "character")
    #env$entrez_genes  = "character",
    #env$uniprot       =  "character",
    #env$go            = "character",
    #env$omim          = "character",
    #env$enzyme        = "character")
)

setGeneric(
  "add<-", 
  function(object, value, ...) StandardGeneric("add<-")
)
setMethod("add<-", c("results"), 
          function(object, value) {
            object@OrderHistory <- append(object@OrderHistory, value)
            object    
          }
)