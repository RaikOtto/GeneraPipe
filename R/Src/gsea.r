eset = eset[ , which( colnames(eset) %in% phenodata$ID  ) ]  
index = which( colnames(phenodata) == cohorts_type )
pData( eset )$Group = phenodata[ match( colnames( eset ), phenodata$ID, nomatch = 0 ), index  ]
#eDatSet = eset

if (var_filter ){
  
  exprs( eDatSet ) = exprs( varFilter( eDatSet  , var.cutoff = 0.6 ) )
}

source("Src/cohort_creation.r")

probeIds = rownames( exprs( eset ) )
binaryGeneSet = rep( 0, length( probeIds ) )
names(binaryGeneSet) = probeIds

mapGeneSet = function( mSigDB_vec, probe_map_hgnc, binaryGeneSet ){
  
  vec = rep("0",length(mSigDB_vec))
  vec[ mSigDB_vec == 1  ] = probe_map_hgnc
  mapping = match( as.character( mSigDB_vec ), probe_map_hgnc, nomatch = 0 )
  res = binaryGeneSet
  res[ which( mapping != 0) ] = 1
  
  return(res)
}

#mapGeneSet = function( mSigDB_vec, geneMap, binaryGeneSet ){

#  mapping = which( mSigDB_vec %in% geneMap$hgnz_symbols )
#  probesOfGeneSet = geneMap$entrez_ids[ mapping ]
#  binaryGeneSet[ probesOfGeneSet ] = 1
#  return( binaryGeneSet )
#}


#mSigDB = read.table("~/Downloads/c2.all.v5.0.symbols.gmt.txt",sep="\t",header=F,fill=T)
#mSigDB = read.table("~/Downloads/h.all.v5.0.symbols.gmt.txt",sep="\t",header=F,fill=T)
#identifier_mSigDB = mSigDB[,c(1,2)]
#mSigDB = mSigDB[,c(-1,-2)]

mSigDB = read.table(c2.all.v5_gsea_file_path, header = F, sep ="\t", fill = T)
#mSigDB = read.table("~/Downloads/c2.cp.v5.0.symbols.gmt.txt",header = F, sep ="\t", fill = T)

mSigDB = mSigDB[,-2]
identifier_mSigDB = rownames(mSigDB)

mSigDB_res = apply( mSigDB[,], FUN = mapGeneSet, hgnc_symbols, binaryGeneSet, MARGIN = 1)

colnames(mSigDB_res) = identifier_mSigDB
mSigDB_res = t(mSigDB_res)

###

library("GSEAlm")
pVals = gsealmPerm( eset, ~Cohort, mSigDB_res, nperm = 100 )

pVals_cor = as.data.frame( apply( pVals, 2, p.adjust, method = "BH", n = nrow(pVals) ) )

pVals_lower = pVals_cor$Lower[ order( pVals_cor$Lower, decreasing = F ) ]
names(pVals_lower) = rownames(pVals_cor)[ order( pVals_cor$Lower, decreasing = F ) ]

pVals_Upper = pVals_cor$Upper[ order( pVals_cor$Upper, decreasing = F ) ]
names(pVals_Upper) = rownames(pVals_cor)[ order( pVals_cor$Upper, decreasing = F ) ]

threshold_lower = quantile( pVals$Lower, probs = seq(0.0,1.0,by = .01))[3]
threshold_upper = quantile( pVals$Upper, probs = seq(0.0,1.0,by = .01))[3]

lower = rownames(pVals)[ pVals$Upper <= threshold_upper ]
higher= rownames(pVals)[ pVals$Lower <= threshold_lower ]

#write.table(pVals, "~/Dropbox/PhD/Ovarian_cancer/Results/pVals.tab",sep="\t",row.names=T,quote=F)
write.xlsx(higher, paste(results_file_path, "upper.xlsx", sep = "/"), row.names = F)
write.xlsx(lower, paste(results_file_path, "lower.xlsx", sep = "/"), row.names = F)
