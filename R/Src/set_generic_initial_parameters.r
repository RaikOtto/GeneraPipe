suppressMessages(library("xlsx"))
suppressMessages(library("gdata"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("limma"))
suppressMessages(library("KEGG.db"))
suppressMessages(library("pathview"))
suppressMessages(library("stringr"))

use_gsea = F
export_eset = F
filter_topall_res = T

output_path = dirname(cel_files_path)
dir.create( paste( output_path, "Output", sep = "/"), showWarnings = F )
cel_files_path = sub( x = cel_files_path, "/$", "")
cpdb_file                   = "CPDB_pathways_genes.tab"
tissue_norm_exprs_file      = "GSE1133-GPL96_series_matrix.txt"
tissue_norm_exprs_file_path = paste( pipeline_loc , paste("Misc", tissue_norm_exprs_file, sep = "/" ), sep = "/")
quality_report_path         = paste( output_path, "Output/QC_report" , sep = "/" )
results_file_path           = paste( output_path, "Output", paste( "Results", project_name, sep ="_"), sep = "/" )
name_res_file               = paste( results_file_path, paste( "dif_exp_results", "csv", sep ="."), sep = "/" )
pathway_maps_path           = paste( results_file_path  , "pathway_maps_dif_exp", sep ="/" )
body_exprs_maps_path        = paste( pipeline_loc, "Misc/HPM_gene_level_epxression_matrix_Kim_et_al_052914.csv", sep ="/" )
absent_gene_file_path       = paste( cel_files_path, absent_genes_file, sep ="/" )
tissue_abbundance_res_file  = paste( cel_files_path, "Tissue_abundance_results.csv", sep ="/" )
kegg_file_path              = paste( cel_files_path , kegg_file, sep ="/")
cpdb_file_path              = paste( pipeline_loc , paste( "Misc" ,cpdb_file, sep ="/" ) , sep ="/" )
c2.all.v5_gsea_file_path    = paste( cel_files_path, "c2.all.v5.1.symbols.gmt", sep = "/")
time_series_res_file_path   = paste( cel_files_path, time_series_res_file, sep ="/")
entities_of_interest_path   = paste( results_file_path, "Entities_of_interest", sep ="/")
genes_of_interest_file_path = paste( entities_of_interest_path, genes_of_interest_file, sep ="/") 
user_folder                 = as.character( system("echo $HOME", intern = T) )
vcf_folder                  = ""
expression_data = "~/Dropbox/PhD/NAR_sub_june_2015/expression_all.txt"
frma_path = "/media/rayott/Backup/Runx_AML1_Leukemea/"

#if (time_series){ quality_control_only = T; qc_control = T }

strEndsWith <- function(haystack, needle)
{
  hl <- nchar(haystack)
  nl <- nchar(needle)
  if(nl>hl)
  {
    return(F)
  } else
  {
    return(substr(haystack, hl-nl+1, hl) == needle)
  }
}

cohorts_file_path = paste( cel_files_path, cohorts_file, sep ="/" )

if (  strEndsWith( cohorts_file_path, ".csv" ) ){  
  phenodata = read.table( cohorts_file_path , header = T , sep = "," )
} else {  
  phenodata = read.table( cohorts_file_path , header = T , sep = "\t" )
}

if (  strEndsWith( kegg_file_path, ".csv" ) ){  
  keggdata = read.table( kegg_file_path , header = T , sep = "," )
} else {  
  keggdata = read.table( kegg_file_path , header = T , sep = "\t" )
}
#cpdbdata = read.table( cpdb_file_path , header = T , sep = "\t", fill =T )
#library("stringr")
#cpdb_id = str_replace(cpdbdata$external_id,"path:","")

#mapping_kegg_cpdb = match(keggdata$HSA_ID, cpdb_id  )
#keggdata$genes = cpdbdata$hgnc_symbol_ids[ mapping_kegg_cpdb  ]

if ( chip_type == "hgu133plus2" ){
  
  suppressMessages(library("hgu133plus2.db"))
  suppressMessages(library("affy"))
  suppressMessages(library("simpleaffy"))
  suppressMessages(library("affyPLM"))
  suppressMessages(library("affycoretools"))
  suppressMessages(library("affyQCReport"))
  suppressMessages(library("annaffy"))
  
} else if ( chip_type == "hgu133a" ){
  
  suppressMessages(library("hgu133a.db"))
  suppressMessages(library("affy"))
  suppressMessages(library("simpleaffy"))
  suppressMessages(library("affyPLM"))
  suppressMessages(library("affycoretools"))
  suppressMessages(library("affyQCReport"))
  suppressMessages(library("annaffy"))
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st", "pd.huex.1.0.st.v2" ) ){
  
  suppressMessages(library("oligoData"))
  suppressMessages(library("pd.huex.1.0.st.v2"))
  
} else if ( chip_type == "HumanHT-12.v4" ){
  
  suppressMessages(library("lumi"))
  suppressMessages(library("GEOquery"))
  
} else {
  
  message(c("Unknown Chip Type: ",chip_type))
  stop()
}

# polishing

cel_files_path = sub( x = cel_files_path, pattern = "~", replacement = u)ser_folder  )
kegg_file_path = sub( x = kegg_file_path, pattern = "~", replacement = user_folder  )
cohorts_file_path = sub( x = cohorts_file_path, pattern = "~", replacement = user_folder  )
quality_report_path = sub( x = quality_report_path, pattern = "//", replacement = "/")
pipeline_loc = sub( x = pipeline_loc, pattern = "//", replacement = "/")
pipeline_loc = sub( x = pipeline_loc, pattern = "~", replacement = "")
