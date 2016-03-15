### Script to invoke limma and GSEA analysis using GeneraPipe. To be called from command line. ###

### Call with: "Rscript GeneraPipeTerminal.R cel_files_path project msigdb_path output_path gsea".
### To add your own project to GeneraPipe database use GeneraPipe/extdata/createDB.R
### For test run call: "Rscript GeneraPipeTerminal.R test GSE29156 test ~/ TRUE".
### Test run will perform analysis on a small test data set. Output will be stored in your home directory, MSigDB will be downloaded from GitHub into ~/Output_GeneraPipe.


# Check if all required packages are installed and can be loaded. If not, install them via Bioconductor or Cran.
packBioc = c("limma", "KEGG.db", "pathview", "hgu133plus2.db", "affy", "affyPLM", "hgu133a.db", "genefilter", "oligo", "biomaRt", "oligoClasses", "Biobase", "pd.huex.1.0.st.v2")
packCran = c("gdata", "stringr", "WriteXLS", "dplyr")

indexBioc = which(! packBioc %in% installed.packages())
indexCran = which(! packCran %in% installed.packages())

if (length(indexBioc) > 0){
  current = packBioc[indexBioc]
  source("https://bioconductor.org/biocLite.R")
  biocLite(current)
}

if (length(indexCran) > 0){
  current = packCran[indexCran]
  install.packages(current)
}

packages = c(packBioc, packCran)

print(suppressMessages(sapply(packages, require, character.only = TRUE, quietly = TRUE)))

library(devtools)

if(!("GeneraPipe" %in% installed.packages())){
  install_github("RaikOtto/GeneraPipe")
}

library(GeneraPipe)

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("At least five arguments must be supplied. 1: cel_files_path, 2: project id, 3: msigdb_path, 4: output_path, 5: gsea TRUE/FALSE", call. = FALSE)
}

cel_files_path = args[1]
project = args[2]
msigdb_path = args[3]
output_path = args[4]
gsea = args[5]

GeneraPipe::StartAnalysis(
  cel_files_path = cel_files_path,
  project = project,
  msigdb_path = msigdb_path,
  output_path = output_path,
  gsea = gsea,
  kegg_for_heatmap = TRUE)