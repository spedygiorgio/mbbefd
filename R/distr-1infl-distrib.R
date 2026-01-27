#d, p, q, r function for one-inflated distribution


doifun <- function(x, dfun, p1, log=FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(x))
  p1 <- rep_len(p1, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  
  #x == 1
  idmain <- x == 1 & p1 >= 0 & p1 <= 1
  if(!log)
  {
    res[idmain] <- p1[idmain]  
  }else
    res[idmain] <- log(p1[idmain])
  # x not in [0,1]
  idmain <- (x > 1 | x < 0) & p1 >= 0 & p1 <= 1
  if(!log)
  {
    res[idmain] <- 0
  }else
    res[idmain] <- -Inf
  # x in [0,1)
  idmain <- 0 <= x & x < 1 & p1 >= 0 & p1 <= 1
  if(!log)
  {
    res[idmain] <- dfun(x[idmain], log=FALSE, ...) * (1 - p1[idmain])
  }else
    res[idmain] <- dfun(x[idmain], log=TRUE, ...) + log(1 - p1[idmain])
  res
}

poifun <- function(q, pfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(q)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(q))
  p1 <- rep_len(p1, length.out=m)
  q <- rep_len(q, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1
  res[idmain] <- pfun(q[idmain], lower.tail = TRUE, log.p = FALSE, ...)*(1 - p1[idmain]) + p1[idmain]*(q[idmain] >= 1)
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}


qoifun <- function(p, qfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(p)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(p))
  p1 <- rep_len(p1, length.out=m)
  p <- rep_len(p, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 < 1 & p >= 0 & p < 1 - p1
  p[idmain] <- p[idmain]/(1-p1[idmain]) #transformed quantile
  if(!lower.tail)
    p[idmain] <- 1-p[idmain]
  if(log.p) 
    p[idmain] <- exp(p[idmain]) 
  
  res[idmain] <- qfun(p[idmain], lower.tail = TRUE, log.p = FALSE, ...)
  
  idmain <- p1 >= 0 & p1 < 1 & p <= 1 & p >= 1 - p1
  res[idmain] <- 1
  res
}  

roifun <- function(n, rfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  if(length(n) > 1)
    n <- length(n)
  if(min(length(p1), n) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), n)
  p1 <- rep_len(p1, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1
  res[idmain] <- rbinom(sum(idmain), 1, p1[idmain])
  
  notone <- res[idmain] != 1 & idmain
  res[notone] <- rfun(sum(notone), ...)
  
  res
}

#exposure curve and moment functions
ecoifun <- function(x, ecfun, mfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(x))
  p1 <- rep_len(p1, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- x >= 0 & x <= 1  & p1 >= 0 & p1 <= 1
  
  G0 <- ecfun(x[idmain], ...) #exposure curve
  E0 <- mfun(order=1, ...) #expectation
  
  res[idmain] <- ((1-p1[idmain])*G0[idmain] + p1[idmain]*x[idmain]/E0)/(1-p1[idmain]+p1[idmain]/E0)
  res
}


# moment function
moifun <- function(order, mfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  stopifnot(is.numeric(order))
  
  if(min(length(p1), length(order)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(order))
  p1 <- rep_len(p1, length.out=m)
  order <- rep_len(order, length.out=m)
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1
  
  E0 <- mfun(order=order[idmain], ...) #expectation
  res[idmain] <- p1[idmain] + (1-p1[idmain])*E0
  res
}

#total loss function
tloifun <- function(p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  if(length(p1) <= 0)
    return(numeric(0))
  
  m <- length(p1)
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1
  res[idmain] <- p1[idmain]
  res
}
