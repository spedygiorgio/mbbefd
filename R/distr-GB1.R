#d, p, q, r function for generalized beta distribution of the first kind
#(no location and scale paramater)


dgbeta <- function(x, shape0, shape1, shape2, log=FALSE)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  
  if(min(length(shape0), length(shape1), length(shape2), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), length(x))
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  id01 <- 0 <= x & x <= 1
  
  if(log)
  {
    idmain <- id01 & shapepos
    res[idmain] <- log(shape0[idmain])+(shape0[idmain]-1)*log(x[idmain])+dbeta(x[idmain]^(shape0[idmain]), 
                                                                               shape1[idmain], shape2[idmain], log=TRUE)
    res[!id01 & shapepos] <- -Inf
  }else
  {
    idmain <- id01 & shapepos
    res[idmain] <- shape0[idmain]*x[id01]^(shape0[idmain]-1)*dbeta(x[idmain]^(shape0[idmain]), 
                                                                   shape1[idmain], shape2[idmain], log=FALSE)
    res[!id01 & shapepos] <- 0
  }
  res
}

pgbeta <- function(q, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  
  if(min(length(shape0), length(shape1), length(shape2), length(q)) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), length(q))
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  q <- rep_len(q, length.out=m)
  
  res <- rep(NaN, m)
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  
  res[shapepos] <- pbeta(q[shapepos]^shape0[shapepos], shape1[shapepos], shape2[shapepos], 
                         lower.tail=TRUE, log.p=FALSE)
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  res
}


qgbeta <- function(p, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  
  if(min(length(shape0), length(shape1), length(shape2), length(p)) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), length(p))
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  p <- rep_len(p, length.out=m)
  
  res <- rep(NaN, m)
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  idmain <- shapepos & 0 <= p & p <= 1
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p)
  res[idmain] <- qbeta(p[idmain], shape1[idmain], shape2[idmain], lower.tail=TRUE, log.p=FALSE)^(1/shape0[idmain])
  res
}  

rgbeta <- function(n, shape0, shape1, shape2)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  if(length(n) > 1)
    n <- length(n)
  
  if(min(length(shape0), length(shape1), length(shape2), n) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), n)
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  res <- rep(NaN, m)
  
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  
  res[shapepos] <- rbeta(sum(shapepos), shape1[shapepos], shape2[shapepos])^(1/shape0[shapepos])
  res
}


ecgbeta <- function(x, shape0, shape1, shape2)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  
  if(min(length(shape0), length(shape1), length(shape2), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), length(x))
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  idmain <- shapepos & 0 <= x & x <= 1
  
  cst2 <- beta(shape1[idmain], 1/shape0[idmain])/beta(shape1[idmain] + shape2[idmain], 1/shape0[idmain])
  res[idmain] <- pbeta(x[idmain]^shape0[idmain], shape1[idmain]+1/shape0[idmain], shape2[idmain])
  res[idmain] <- res[idmain] + x[idmain]*(1 - pbeta(x[idmain]^shape0[idmain], shape1[idmain], shape2[idmain]))*cst2
  res
}

mgbeta <- function(order, shape0, shape1, shape2)
{
  #sanity check
  stopifnot(is.numeric(shape0))
  stopifnot(is.numeric(shape1))
  stopifnot(is.numeric(shape2))
  
  if(min(length(shape0), length(shape1), length(shape2), length(order)) <= 0)
    return(numeric(0))
  
  m <- max(length(shape0), length(shape1), length(shape2), length(order))
  shape0 <- rep_len(shape0, length.out=m)
  shape1 <- rep_len(shape1, length.out=m)
  shape2 <- rep_len(shape2, length.out=m)
  order <- rep_len(order, length.out=m)
  
  res <- rep(NaN, m)
  
  shapepos <- shape0 > 0 & shape1 > 0 & shape2 > 0
  
  idmain <- shapepos & 0 < shape1 + order / shape0
  res[idmain] <- beta(shape1[idmain] + order[idmain] / shape0[idmain], shape2[idmain]) / beta(shape1[idmain], shape2[idmain])
  res
}
  

###################
#internal functions

#incomplete beta function
betainc <- function(x, a, b) pbeta(x, a, b) * beta(a, b)


#Theil index, see package ineq for other income index (e.g. Gini coefficient)
Theil.theo  <- function(shape0, shape1, shape2)
{
  EX <- beta(shape1+shape2, 1/shape0) / beta(shape1, 1/shape0)
  1/shape0*(digamma(shape1+1/shape0) - digamma(shape1+shape2+ 1/shape0)) - log(EX)
}

Theil.theo.shape0  <- function(shape0, obs)
{
  #compute shape1/shape2 on a rescaled sample and moment estimator
  obs <- obs^shape0
  n <- length(obs)
  m <- mean(obs)
  v <- (n - 1)/n*var(obs)
  aux <- m*(1-m)/v - 1
  shape1 <- m*aux
  shape2 <- (1-m)*aux
  
  Theil.theo(shape0, shape1, shape2)
}
