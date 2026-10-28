# File overview: random generation of the MBBEFD distribution.
# For a single parameter value the fast C++ generators of src/mbbefdRcppCode.cpp
# are used; when the parameters are vectors (one value per variate), the R
# versions rmbbefdR() / rMBBEFDR() of distr-mbbefdR-*.R are used instead.
# The d, p, q, ec, m, tl functions are in distr-mbbefdR-1stparam.R (a,b) and
# distr-mbbefdR-2ndparam.R (g,b).
#' @useDynLib mbbefd
#' @importFrom Rcpp sourceCpp

#see d, p, q, tl, m functions in distr-mbbefdR*.R

#random generation
# internal wrappers around the C++ functions registered by Rcpp (see RcppExports.R)
rmbbefdCpp <- function(n, a, b)
{
  .rmbbefdC(n,  a,  b)
}

rMBBEFDCpp <- function(n, g, b)
{
  .rMBBEFDC(n, g, b)
}


### r function MBBEFD(a,b) for users
# User-level generator for MBBEFD(a,b): n variates on [0,1].
rmbbefd <- function(n, a, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  # a single n must be a non-negative integer (n > 2^31 - 1 cannot be passed to
  # the C++ generators), as in the r*() functions of base R
  if (length(n) == 1 && (is.na(n) || n < 0 || n > .Machine$integer.max))
    stop("invalid arguments")

  # an empty argument gives an empty result
  if (min(length(a), length(b), length(n)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(n))
  # scalar parameters: C++ generator; otherwise: vectorized R generator
  if (m == 1)
    res <- rmbbefdCpp(n, a, b)
  else
    res <- rmbbefdR(n, a, b)
  res
}

### r function MBBEFD(g,b) for users
# User-level generator for MBBEFD(g,b): n variates on [0,1].
rMBBEFD <- function(n, g, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))
  # a single n must be a non-negative integer (n > 2^31 - 1 cannot be passed to
  # the C++ generators), as in the r*() functions of base R
  if (length(n) == 1 && (is.na(n) || n < 0 || n > .Machine$integer.max))
    stop("invalid arguments")

  if (min(length(g), length(b), length(n)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(n))

  # scalar parameters: C++ generator; otherwise: vectorized R generator
  if (m == 1)
    res <- rMBBEFDCpp(n, g, b)
  else
    res <- rMBBEFDR(n, g, b)
  res
}
