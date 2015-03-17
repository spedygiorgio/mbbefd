

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
    res[x < 0] <- 0
    res[x > 1] <- 1
  }else
  {
    res <- log((a+b^x)/(a+1))/log((a+b)/(a+1))  
    res[x < 0] <- 0
    res[x > 1] <- 1  
  }
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
	
	
### R version of d,p,q,r functions MBBEFD(g,b)


		

dMBBEFDR <- function(x, g, b, log=FALSE)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(x)))
  
  if(g == 1 || b == 0) #Dirac
  {
    res <- rep(0, length(x))
    res[x == 1] <- 1
  }else if(g == 1/b && b < 1) #bg=1
  {
    res <- -log(b)*b^x
    res[x == 1] <- b
  }else if(g > 1 && b == 1) #b=1
  {
    res <- (g-1)/(1+(g-1)*x)^2
    res[x == 1] <- 1/g
  }else
  {
    res <- -(1-b)*log(b)*b^(1-x)/((g-1)*b^(1-x) + 1-g*b)^2
    res[x == 1] <- 1/g
  }
  res[x > 1] <- 0
  res[x < 0] <- 0
  
  if(log)
    res <- log(res)
  return(res)  
}  

pMBBEFDR <- function(q, g, b, lower.tail = TRUE, log.p = FALSE)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(q)))
  
  if(g == 1 || b == 0) #Dirac
  {
    res <- rep(0, length(q))
  }else if(g == 1/b && b < 1) #bg=1
  {
    res <- 1-b^q
  }else if(g > 1 && b == 1) #b=1
  {
    res <- 1-1/(1+(g-1)*q)
  }else
  {
    res <- 1- (1-b)/((g-1)*b^q+1-g*b)
  }
  res[q >= 1] <- 1
  res[q < 0] <- 0
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  return(res)  
}  


qMBBEFDR <- function(p, g, b, lower.tail = TRUE, log.p = FALSE)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(p)))
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  if(g == 1 || b == 0) #Dirac
  {
    res <- rep(0, length(p))
    res[p > 0] <- 1
  }else if(g == 1/b && b < 1) #bg=1
  {
    res <- log(1-p)/log(b)
    res[p > 1-b] <- 1
  }else if(g > 1 && b == 1) #b=1
  {
    res <- p/((1-p)*(g-1))
    res[p > 1-1/g] <- 1
  }else
  {
    res <- 1- log((g*b-1)/(g-1) + (1-b)/((1-p)*(g-1)))/log(b)
    res[p > 1-1/g] <- 1
  }
  res[p < 0 | p > 1] <- NaN
  
}  


rMBBEFDR <- function(n, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, n))
  
}


gMBBEFDR <- function(x, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(x)))
  
}

mMBBEFDR <- function(order, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(order)))
  
}
