library("pathview")
#print("Step 10: Creating Pathway maps")

if ( !dir.exists( pathway_maps_path )   ){ dir.create(pathway_maps_path) }
if ( ! exists("eset")  )
  source("Src/normalization.r")

split_fun = function( entry,pos ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[pos] ) } else{ return( "" ) } }

if ( chip_type %in% c( "hgu133plus2" ) ){
  
  library("hgu133plus2.db")
  hgnc_symbols = as.character( unlist( mget( rownames(volc_all), hgu133plus2SYMBOL) ))

} else if ( chip_type %in% c( "hgu133a" ) ){
  
  library("hgu133a.db")
  hgnc_symbols = as.character( unlist( mget( rownames(volc_all), hgu133aSYMBOL) ) )
  
} else if ( chip_type %in% c( "HumanHT-12.v4") ){
  
  hgnc_symbols = ncbifd$Gene.symbol
  hgnc_symbols = hgnc_symbols[ 93:length( hgnc_symbols ) ]

} else if ( chip_type %in% c( "pd.hugene.2.0.st", "pd.huex.1.0.st.v2" ) ){
  
  if ( ( ! ( "HGNC_symb" %in% colnames( topall_res ) ) )  |  ( length( topall_res$HGNC_symb ) == 0 )  ){
    
    print("Pathway map creation error: topall_res$HGNC_symb does no exist")
    quit()
  }

  hgnc_symbols = str_trim( unlist( lapply( volc_all$geneassignment, FUN=split_fun,2 ) ) )
  hgnc_names   = str_trim( unlist( lapply( volc_all$geneassignment, FUN=split_fun,3 ) ) )
} else{ 

  # note: check if this works for affy 
  hgnc_symbols = str_trim( unlist( lapply( fData(eset)$geneassignment, FUN=split_fun,2 ) ) )
  hgnc_names   = str_trim( unlist( lapply( fData(eset)$geneassignment, FUN=split_fun,3 ) ) )  
}

library("biomaRt")

ensembl     = useMart("ensembl",dataset="hsapiens_gene_ensembl")
entrez_ids  = getBM( attributes = c( "entrezgene", "hgnc_symbol" ), values = unique( hgnc_symbols), filters = "hgnc_symbol" , mart = ensembl)

hgnc_map    = match( hgnc_symbols, entrez_ids$hgnc_symbol  , nomatch = 0 )
entrez      = rep( "", length( hgnc_symbols ) )
entrez[ which( hgnc_map != 0 ) ] = entrez_ids$entrezgene[ hgnc_map ]
entrez[ is.na(entrez) ] = ""

if ( chip_type == "HumanHT-12.v4" ){
  
  probe_ids = ncbifd$ID
  eset_select = exprs(eDatSet)[ match( probe_ids, rownames( exprs( eset ) )  ) , ]
  eset_select = eset_select[ which( !is.na( rownames(eset_select) ) ), ]
  exprs_case = round( rowMeans( eset_select[ ,index_case ] ), 2 )
  exprs_ctrl = round( rowMeans( eset_select[ ,index_ctrl ] ), 2 )
  dif_exp = round( exprs_case - exprs_ctrl, 2)
  
} else{
  
  exprs_case = round( rowMeans( exprs( eDatSet )[ ,index_case ] ), 2 )
  exprs_ctrl = round( rowMeans( exprs( eDatSet )[ ,index_ctrl ] ), 2 )
  dif_exp = round( exprs_case - exprs_ctrl, 2)
}

res_all_path = data.frame(

  "logFC"               = dif_exp,
  "HGNC"                = hgnc_symbols,
  "entrez"              = entrez
)

setwd( pathway_maps_path )
Kegg_id = as.character( keggdata$Kegg_id )

ent = res_all_path$entrez
ent = ent[ent != ""]
exp_dat = res_all_path$logFC[  res_all_path$entrez != "" ]
names(exp_dat) = res_all_path$entrez[  res_all_path$entrez != "" ]

kegg <- mget( as.character( ent  ), KEGGEXTID2PATHID, ifnotfound=list(NA))

for (id in Kegg_id  ){
  
  pathview( gene.data = exp_dat, pathway.id = id, kegg.native = T, limit = max( abs(exp_dat)  ) )
}

system("rm *.xml")
system("ls | grep -v 'pathview.png$' | xargs rm")

setwd( pipeline_loc )
