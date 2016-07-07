### Script in development
### Use this script to add your own project to the GeneraPipe database. GeneraPipe needs to be installed("install_github("RaikOtto/GeneraPipe")").
### Add all the parameters to the data.frame "newParameter" and execute the script.
### To run GeneraPipe with your new project, call StartAnalysis(cel_files_path, project = Your_ID, msigdb_path, output_path, gsea)

library(dplyr)
package_path = system.file("", package = "GeneraPipe")
setwd(package_path)
defaultDB = src_sqlite("GeneraPipeDefaultDB.sqlite3")

newParameters = data.frame( 
  ID           = c("GSE21344"), # ID to access your project. Example: "GSE29156", needs to be of class character
  project_name = c("GSE21344_goldenSpike"), # Project name. Example: "GSE29156_test", needs to be of class character
  chip_type    = c("drosophila2"), # chip type used for experiment. Example: "pd.huex.1.0.st.v2", needs to be of class character
  cohorts_file = c("cohorts.tab"), # name of cohort file. Example: "cohorts.tab", needs to be of class character
  set_ctrl     = c("condition_A"), # identifier for control cohort in cohort-file. Example: "Normal", needs to be of class character
  set_case     = c("condition_B"), # identifier for control cohort in cohort-file. Example: "Tumor", needs to be of class character
  kegg_file    = c("kegg_pathways_of_interest.csv"), # name of kegg pathway file. Example: "pathways.csv", needs to be of class character
  zipped       = c(TRUE),   # are the CEL files zipped? Example: TRUE, needs to be of class logical
  p_val        = c(0.05),   # p value threshold for diff expression analysis. Example: 0.05, needs to be of class integer
  lfc_exp      = c(1)   # log fold change threshold for diff expression analysis. Example: 1, needs to be of class integer
)

db_insert_into( defaultDB$con, table = "projects", values = newParameters)
