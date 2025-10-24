#' Binary to decimal helper function
#'
#' @param x string binary value to convert
#' @return numeric decimal value
#' @examples
#' BinToDec("1001")
#' @export
BinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

#' Decimal to binary helper function
#'
#' @param x numeric value to convert
#' @param len numeric length of required binary string
#' @return string binary value
#' @examples
#' DecToBin(9, 4)
#' @export
DecToBin <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(paste(s, collapse=""))
}

#' Decimal to (numeric) binary helper function
#'
#' @param x numeric value to convert
#' @param len numeric length of required binary string
#' @return numeric binary vector
#' @examples
#' DecToBinV(9, 4)
#' @export
DecToBinV <- function(x, len) {
  s = c()
  for(j in (len-1):0)
  {
    if(x >= 2**j) { s=c(s,1); x = x-2**j } else { s=c(s,0)}
  }
  return(s)
}


#' Difference helper function for vectorised calculations
#'
#' @param x numeric
#' @param y numeric
#' @return numeric difference
#' @examples
#' matdiff(5, 2)
#' @export
matdiff <- function(x, y) {
  return(x-y)
}

#' Check compatibility between two bitstrings, reported as a distance
#'
#' @param s1 character bit string 1
#' @param s2 character bit string 2
#' @return numeric difference
#' @examples
#' is_compat("101", "001")
#' @export
is_compat = function(s1, s2) {
  diffs = as.numeric(strsplit(s2, "")[[1]]) - as.numeric(strsplit(s1, "")[[1]])
  if(min(diffs) < 0) {
    return(Inf)
  } else if(sum(diffs) == 0) {
    return(Inf)
  } else {
    return(sum(diffs))
  }
}

#' Check compatibility between two numeric bitstrings, reported as a distance
#'
#' @param s1 numeric bit vector 1
#' @param s2 numeric bit vector 2
#' @return numeric difference
#' @examples
#' is_compat_num(c(1,0,1), c(0,0,1))
#' @export
is_compat_num = function(s1, s2) {
  diffs = s2-s1
  sumdiffs = sum(diffs)
  if(min(diffs) < 0) {
    return(Inf)
  } else if(sumdiffs == 0) {
    return(Inf)
  } else {
    return(sumdiffs)
  }
}
