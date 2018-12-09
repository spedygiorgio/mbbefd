

### R version of d,p,q,r functions MBBEFD(a,b)

dmbbefdR <- function(x, a, b, log=FALSE)
{
  if(is.infinite(b))
    return(rep(NaN, length(x)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(x)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(x)))
  }
  
  if(a == 0 || b == 1 || a == -1) #Dirac
  {
    res <- rep(0, length(x))
    res[x == 1] <- 1
  }else if(is.infinite(a))
  {
    res <- -log(b)*b^x
  }else
  {
    res <- -a * (a+1) * b^x * log(b) / (a + b^x)^2 
    res[x == 1] <- (a+1) * b / (a+b)
  }
  res[x > 1] <- 0
  res[x < 0] <- 0
  
  if(log)
    res <- log(res)
  res  
}  
	
pmbbefdR <- function(q, a, b, lower.tail = TRUE, log.p = FALSE)
{
  
  if(is.infinite(b))
    return(rep(NaN, length(q)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(q)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(q)))
  }
  
  if(a == 0 || b == 1 || a == -1) #Dirac
  {
    res <- rep(0, length(q))
  }else if(is.infinite(a))
  {
    res <- 1-b^q
  }else
  {
    res <- a * ( (a+1) / (a + b^q) - 1) 
  }
  res[q >= 1] <- 1
  res[q <= 0] <- 0 #modified 9/12/2018

  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}  
  
	 
qmbbefdR <- function(p, a, b, lower.tail = TRUE, log.p = FALSE)
{
  if(is.infinite(b))
    return(rep(NaN, length(p)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(p)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(p)))
  }
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  if(a == 0 || b == 1 || a == -1) #Dirac
  {
    res <- rep(0, length(p))
    res[p > 0] <- 1
  }else if(is.infinite(a))
  {
    res <- rep(1, length(p))
    res[p < 1-b] <- log(1-p[p < 1-b])/log(b)
  }else
  {
    pab <- (a+1)*b/(a+b)
    res <- rep(1, length(p))
    p2 <- p[p < 1-pab]
    res[p < 1-pab] <- log((1-p2)*a/(a+p2))/log(b)
  }
  res[p < 0 | p > 1] <- NaN
  
  res
}  

  
rmbbefdR <- function(n, a, b)
{
  if(is.infinite(b))
    return(rep(NaN, length(n)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(n)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(n)))
  }

  qmbbefdR(runif(n, 0, 1), a, b)
}
	
	
ecmbbefdR <- function(x, a, b)
{
  if(is.infinite(b))
    return(rep(NaN, length(x)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(x)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(x)))
  }
  
  if(a == 0 || b == 1 || a == -1) #Dirac
  {
    res <- x
  }else if(is.infinite(a))
  {
    res <- (1-b^x)/(1-b)
  }else
  {
    res <- log((a+b^x)/(a+1))/log((a+b)/(a+1))  
  }
  res[x < 0] <- 0
  res[x > 1] <- 1
  
  res
}
	
#moment
mmbbefdR <- function(order, a, b)
{
  if(is.infinite(b))
    return(rep(NaN, length(order)))
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(rep(NaN, length(order)))
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(rep(NaN, length(order)))
  }
  if(any(order > 2))
    stop("not yet implemented.")
  
  if(a == 0 || b == 1 || a == -1) #Dirac
  {
    res <- rep(1, 2)
  }else if(is.infinite(a))
  {
    res <- c((b-1)/log(b), 2*pgamma(-log(b),2)*gamma(2)/log(b)^2)
  }else
  {
    res <- c(log((a+b)/(a+1))/log(b)*(a+1), 2*(a+1)/log(b)*(log(a+b) - gendilog(a,b)))
  }
  return(res[order])
}
	
#total loss
tlmbbefdR <- function(a, b)
{
  if(is.infinite(b))
    return(NaN)
  if(is.infinite(a))
  {  
    if(!(b > 0 && b < 1))
      return(NaN)
  }else if(a != -1)
  {
    if(!(a +1 >0 && b > 0 && a*(1-b) >= 0))
      return(NaN)
  }
  
  if(is.infinite(a))
    res <- b
  else
    res <- (a+1)*b/(a+b)
  res
}



### d,p,q,ec,m,tl functions MBBEFD(a,b)

dmbbefd <- dmbbefdR

pmbbefd <- pmbbefdR

qmbbefd <- qmbbefdR

ecmbbefd <- ecmbbefdR

mmbbefd <- mmbbefdR

tlmbbefd <- tlmbbefdR

