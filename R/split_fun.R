split_fun = function( entry, pos ){ 
  res = unlist( stringr::str_split( entry, " // " ) )
  if (length(res) > 1){
    return( res[pos] ) 
  } else{ 
    return( "" ) 
  }
}