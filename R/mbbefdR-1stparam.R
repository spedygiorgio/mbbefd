

### R version of d,p,q,r functions MBBEFD(a,b)

dmbbefdR <- function(x, a, b, log=FALSE)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(x)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(x))
    res[x == 1] <- 1
  }else
  {
    res <- -a * (a+1) * b^x * log(b) / (a + b^x)^2 
    res[x == 1] <- (a+1) * b / (a+b)
  }
  res[x > 1] <- 0
  res[x < 0] <- 0
  
  if(log)
    res <- log(res)
  return(res)  
}  
	
pmbbefdR <- function(q, a, b, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(q)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(q))
    res[q == 1] <- 1
  }else
  {
    res <- a * ( (a+1) / (a + b^q) - 1) 
    res[q >= 1] <- (a+1) * b / (a+b)
    res[q < 0] <- 0
  }
  res[q >= 1] <- 1
  res[q < 0] <- 0

  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  return(res)  
}  
  
	 
qmbbefdR <- function(p, a, b, lower.tail = TRUE, log.p = FALSE)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(p)))
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(p))
    res[p > 0] <- 1
  }else
  {
    pab <- (a+1)*b/(a+b)
    res <- log((1-p)*a/(a+p))/log(b)
    res[p >= 1-pab] <- 1
  }
  res[p < 0 | p > 1] <- NaN
  
  return(res)  
}  

  
rmbbefdR <- function(n, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, n))
  qmbbefdR(runif(n, 0, 1), a, b)
}
	
	
gmbbefdR <- function(x, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(x)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- x
  }else
  {
    res <- log((a+b^x)/(a+1))/log((a+b)/(a+1))  
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}
	
mmbbefdR <- function(order, a, b)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
    return(rep(NaN, length(order)))
  
  if(order == 1)
    return(log((a+b)/(a+1))/log(b)*(a+1))
  else
    stop("not yet implemented.")
}
	
	