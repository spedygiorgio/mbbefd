#d, p, q, r function for generalized beta distribution of the first kind
#(no location and scale paramater)

#should it be dstpareto01?
dgbeta <- function(x, shape0, shape1, shape2, log=FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(x)))
  
  res <- 1/shape0*x^(1/shape0-1)*dbeta(x^(1/shape0), shape1, shape2, log=log)
  res[x == 0] <- 0
  res
}

pgbeta <- function(q, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(q)))
  pbeta(q^(1/shape0), shape1, shape2, lower.tail=lower.tail, log.p=log.p)
}


qgbeta <- function(p, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(p)))
  qbeta(p, shape1, shape2, lower.tail=lower.tail, log.p=log.p)^(shape0)
}  

rgbeta <- function(n, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, ifelse(length(n)>1, length(n), n)))
  rbeta(n, shape1, shape2)^(shape0)
}

#internal function
betainc <- function(x, a,b) pbeta(x, a, b)*beta(a,b)

ecgbeta <- function(x, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(x)))
  
  b12 <- beta(shape1, shape2)
  EX <- mgbeta(1, shape0, shape1, shape2)
  
  betainc(x^(1/shape0), shape1+shape0, shape2)/(b12*EX) + x*(b12 - betainc(x^(1/shape0), shape1, shape2))/(b12*EX)
}

mgbeta <- function(order, shape0, shape1, shape2)
{
  if(!is.numeric(shape0) || !is.numeric(shape1) || !is.numeric(shape0))
    stop("non numeric argument.")
  if(shape0 < 0 || shape1 < 0 || shape2 < 0)
    return(rep(NaN, length(order)))
  
  beta(shape1+shape2, order*shape0) / beta(shape1, order*shape0)
}
  
