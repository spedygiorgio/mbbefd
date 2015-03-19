#d, p, q, r function for one-inflated shifted truncated Pareto distribution


doistpareto <- function(x, a, p1, log=FALSE)
{
  if(!(a > 0 && p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(x)))
  
  res <- dstpareto(x, a, log=FALSE)*(1 - p1)
  res[x == 1] <- p1
  
  if(log)
    res <- log(res) 
  res
}

poistpareto <- function(q, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a > 0 && p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(q)))
  
  res <- pstpareto(q, a, lower.tail = TRUE, log.p = FALSE)*(1 - p1)
  res[x >= 1] <- 1
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}


qoistpareto <- function(p, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a > 0 && p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(p)))
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  res <- qstpareto(p/(1-p1), a, lower.tail = TRUE, log.p = FALSE)
  res[p >= 1-p1] <- 1
  
  res
}  

roistpareto <- function(n, a, p1)
{
  if(!(a > 0 && p1 >= 0 && p1 <= 1))
    return(rep(NaN, n))
  qoistpareto(runif(n, 0, 1), a, p1)
}


eoicstpareto <- function(x, a, p1)
{
  if(!(a > 0 && p1 >= 0 && p1 <= 1))
    return(rep(NaN, length(x)))
  
  if(a == 1)
  {
    res <- (2*log(x+1) - x)/(2*log(2) - 1)
  }else
  {
    res <- ((x+1)^(-a+1) - 2^(-a)*x*(-a+1) - 1)/(2^(-a+1)-2^(-a)*(-a+1) - 1)
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}
