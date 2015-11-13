#print("Comparison to normal tissue expression")

library("hgu133a.db")
mapWithMultiProbes_symbols = toggleProbes( hgu133aSYMBOL, "all")

if(F){
  tissue_of_interest = c("kidney","kidney.1")
  
  data_tmp = read.table(tissue_norm_exprs_file_path, sep ="\t",header=T)
  tissue_norm_exprs = data_tmp[ , 2 : length( data_tmp[ 1 , ] ) ] ; remove( data_tmp )
  tissue_norm_exprs = log2(tissue_norm_exprs)
  
  print(c("Position of samples of interest:", which( colnames(tissue_norm_exprs) %in% tissue_of_interest  ) ))
  
  probe_ids = read.table(tissue_norm_exprs_file_path, sep ="\t",header=T)[,1]
  tissue_hgnc_symbs = mget( as.character( probe_ids ), mapWithMultiProbes_symbols )
  
  print(c("Amount probes contained in both topall and tissue expression data:",sum(rownames(topall) %in% probe_ids)))
}
#tmp = as.matrix( tissue_norm_exprs, ncol = 158 )
#colnames( tmp ) = colnames(tissue_norm_exprs)
#tissue_norm_exprs = tmp; remove(tmp)

topall_high = topall[topall$logFC > 0,]

if ( chip_type %in% c("hgu133plus2" ) ){
  
  hgnc_high_topall = mget( rownames(topall_high), hgu133plus2SYMBOL )

} else if ( chip_type %in% c("hgu133a", "pd.hg.u133a" ) ){
  
  library("hgu133a.db")
  hgnc_high_topall = mget( rownames(topall_high), hgu133aSYMBOL )

} else if ( chip_type %in% c( "pd.hugene.2.0.st", "pd.huex.1.0.st.v2" ) ){
  
  featureData( eset  ) = getNetAffx( eset, type = "transcript" )
  split_fun = function( entry ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[2] ) } else{ return( "" ) } }
  
  mapping = match( rownames(topall_high), rownames(eset )  )
  sub_set = eset[ mapping ,]
  hgnc_high_topall = str_trim( unlist( lapply( featureData( sub_set  )$geneassignment, FUN=split_fun ) ) )
  hgnc_genes = str_trim( unlist( lapply( featureData( eset  )$geneassignment, FUN=split_fun ) ) )
  remove(sub_set)
  
}

#mapping = match( hgnc_topall, tissue_hgnc_symbs, nomatch = 0  )

#subset = tissue_norm_exprs[ mapping ,]
#res = data.frame( 
#  "hgnc_topall" = as.character( tissue_hgnc_symbs[ mapping ] ), 
#  "row_means_value" = rowMeans(subset) 
#)

#res = res[ res$hgnc_symbol != "NA" ,]
#res = res[order( as.double( res$row_means_value), decreasing = F ),]

### run blood 

print("Calculating blood genes")

if ( !exists("eset_blood")  ){
  
  library("pd.hg.u133a")
  celFiles    = list.celfiles( "~/Dropbox/PhD/Kidney_Cancer_Biomarker/GSE2888_RAW_GPL96", full =T, listGzipped = T )
  raw_data_blood= read.celfiles( filenames = celFiles) # type = pd.hg.u133a
  eset_blood  = rma( raw_data_blood, normalize = T)
  blood_genes = mget( rownames( eset_blood ), hgu133aSYMBOL )
}

### proteome map

expr_map = read.table( body_exprs_maps_path, sep =",",header=T)
expr_map = expr_map[ , - which( startsWith( colnames(expr_map), "Fetal" ) ) ]

expr_gene  = expr_map[, 1]
expr_map   = expr_map[,-1]

### merge results

eset_mapping      = match( rownames( topall_high ), rownames( eset ), nomatch = 0)
mapping_eset      = which( rownames( topall_high ) %in% rownames( eset ) )
hgnc_high_topall  = hgnc_high_topall[ !is.na(hgnc_high_topall)  ]

topall_high_filt = topall_high[ which( as.character( hgnc_high_topall ) %in% as.character( expr_gene)  )  ,]

mapping_blood = match( rownames( topall_high_filt ), rownames( eset_blood ), nomatch = 0 )
topall_high_filt = topall_high_filt[ which( rownames(topall_high_filt) %in%  rownames( eset_blood ) ),]
eset_tmp = eset[ which( rownames(eset) %in% rownames( eset_blood ) )  ,]
hgnc_genes_filt = hgnc_genes[ match( rownames(topall_high_filt), rownames( eset_tmp ), nomatch = 0) ]

topall_high_filtmap = match( rownames(topall_high_filt), rownames(eset_tmp), nomatch = 0  )

res = cbind(
  
  hgnc_genes_filt,
  round( topall_high_filt$logFC ,2 ),
  round( rowMeans(
    exprs( eset_tmp[ topall_high_filtmap, ] )
  ), 2),
  round( rowMeans(
    exprs( eset_blood[ mapping_blood, ] )
  ), 2)
)

remove(eset_tmp)

expr_map_mapping = match( hgnc_genes_filt , as.character( expr_gene ), nomatch = 0 )
expr_map_present = which( hgnc_genes_filt %in% as.character( expr_gene ) )

res = res[ expr_map_present  ,]
tmp = apply( expr_map[ expr_map_mapping, ], MARGIN = 2, FUN = function(x){ x[ x < 1.0] = 1.0; return(  x )  } )
tmp = round( log2( tmp  ), 1)

counter = apply( tmp, MARGIN= 1, FUN = function(x){ return( sum( x >= 4  ) ) }   )

res = cbind( res, counter, tmp  )

colnames(res)[1:5] = c("HGNC_symbol","LogFC","Exprs_Strength", "Exprs_Blood","NMBR_Exprs_Tissue_Greater_4")
res = res[ order( unlist( res[,5]), decreasing = F ) , ]

write.table( res, tissue_abbundance_res_file , row.names=F, quote = F, sep =",", col.names =T  )

l = length( res[1,] )
my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = l)
rownames(tmp) = res[,1]

library("gplots")

png(paste(cel_files_path,"tissue_abbundace.png", sep ="/"))
  heatmap.2(t(tmp), col=my_palette, trace="none")
dev.off()
