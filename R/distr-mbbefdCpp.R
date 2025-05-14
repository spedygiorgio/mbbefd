#' @useDynLib mbbefd
#' @importFrom Rcpp sourceCpp

#see d, p, q, tl, m functions in distr-mbbefdR*.R

#random generation
rmbbefdCpp <- function(n, a, b)
{
  .rmbbefdC(n,  a,  b) 
}

rMBBEFDCpp <- function(n, g, b)
{
  .rMBBEFDC(n, g, b) 
}
  

### r function MBBEFD(a,b) for users
rmbbefd <- function(n, a, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(n)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(n))
  if(m == 1)
    res <- rmbbefdCpp(n, a, b)
  else
    res <- rmbbefdR(n, a, b)
  res
}

### r function MBBEFD(g,b) for users
rMBBEFD <- function(n, g, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))
  
  if(min(length(g), length(b), length(n)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(n))
  
  if(m == 1)
    res <- rMBBEFDCpp(n, g, b)
  else
    res <- rMBBEFDR(n, g, b)
  res
}

