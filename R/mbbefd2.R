

### check 

dMBBEFD <- function(x, a, b)
{
  if(a +1 >0 && b > 0 && a*(1-b) >= 0)
    return(rep(NaN, length(x)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(x))
    res[x == 1] <- 1
    return(res)
  }
  
  res <- -a * (a+1) * b^x * log(b) / (a + b^x)^2 
  res[x == 1] <- (a+1) * b / (a+b)
  return(res)  
}  
	
pMBBEFD <- function(q, a, b)
{
  if(a +1 >0 && b > 0 && a*(1-b) >= 0)
    return(rep(NaN, length(q)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(q))
    res[q == 1] <- 1
    return(res)
  }
  
  res <- a * ( (a+1) / (a + b^q) - 1) 
  res[q >= 1] <- (a+1) * b / (a+b)
  res[q < 0] <- 0
  return(res)  
}  
  
	 
qMBBEFD <- function(p, a, b)
{
  if(a +1 >0 && b > 0 && a*(1-b) >= 0)
    return(rep(NaN, length(p)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- rep(0, length(p))
    res[p > 0] <- 1
    return(res)
  }
  
  pab <- (a+1)*b/(a+b)
  
  res <- log((1-p)*a/(a+p))/log(b)
  res[p >= 1-pab] <- 1
  res[p < 0 | p > 1] <- NaN
  return(res)  
}  

  
	
	
gMBBEFD <- function(x, a, b)
{
  if(a +1 >0 && b > 0 && a*(1-b) >= 0)
    return(rep(NaN, length(x)))
  
  if(a == 0 || b == 1) #Dirac
  {
    res <- x
    res[x < 0] <- 0
    res[x > 1] <- 1
    return(res)
  }
  
  res <- log((a+b^x)/(a+1))/log((a+b)/(a+1))	
  res[x < 0] <- 0
  res[x > 1] <- 1
  res
}
	
mMBBEFD <- function(a,b)
	log((a+b)/(a+1))/log(b)*(a+1)
	
		