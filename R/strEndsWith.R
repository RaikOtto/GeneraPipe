strEndsWith <- function(haystack, needle)
{
  hl <- nchar(haystack)
  nl <- nchar(needle)
  if(nl > hl)
  {
    return(FALSE)
  } else
  {
    return(substr(haystack, hl-nl+1, hl) == needle)
  }
}