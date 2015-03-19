

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
  res
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
  res
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
  
  res
}  


rMBBEFDR <- function(n, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, n))
  qMBBEFDR(runif(n, 0, 1), g, b) 
}


ecMBBEFDR <- function(x, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(x)))
  
  if(g == 1 || b == 0) #Dirac
  {
    res <- x
  }else if(g == 1/b && b < 1) #bg=1
  {
    res <- (1-b^x)/(1-b)
  }else if(g > 1 && b == 1) #b=1
  {
    res <- log(1+(g-1)*x)/log(g)
  }else
  {
    res <- log((g-1)*b/(1-b)+(1-g*b)*b^x/(1-b))/log(g*b)
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}

#moment
mMBBEFDR <- function(order, g, b)
{
  if(!(g >= 1 && b >= 0))
    return(rep(NaN, length(order)))
  0
}


#total loss
tlMBBEFDR <- function(g, b)
{
  if(!(g >= 1 && b >= 0))
    return(NaN)
  1/g
}

