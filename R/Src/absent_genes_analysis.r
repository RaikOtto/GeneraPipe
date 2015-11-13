suppressMessages(library("genefilter"))

f1 = pOverA( 1/3, log2( 100 ) )
f2 = function( x ) ( diff( range( x, na.rm = T ) ) > log2( 1.5 ) )
ff = filterfun( f1, f2 )

not_absent_index = genefilter( eset, ff )
eDatSet = eset[ not_absent_index, ]  #filtered expression matrix

probes_included = rownames( eset[ not_absent_index, ] )
probes_excluded = rownames( eset[ !not_absent_index, ] )

gene_list = unlist( str_split( as.character(keggdata$Gene_id_hgnc ), "," ) )
#gene_list = unique( c( genes_of_interest,gene_list ) )

genes_included = gene_list[ not_absent_index ]
genes_excluded = gene_list[ !not_absent_index ]

#mapping_interesting = which( genes_excluded %in% unlist(str_split( interesting_pathways_table$hgnc_symbol_ids,",") ) )

#absent_report = data.frame(

 # "probe_id" = probes_excluded[ mapping_interesting ],
#  "hgnc_gene" = genes_included,
 # "expression" = round( exprs(eset)[ mapping_interesting, ], 2 )
#)

if (exists("absent_report")){
  colnames( absent_report ) = c("probe_id","hgnc_gene",colnames(exprs(eset)[ - not_absent_index, ]))
  write.table( absent_report, absent_gene_file_path, row.names=F, sep = "\t", quote = F )
  remove(absent_report)
}