
fitDR <- function(x, dist, method="mle", ...)
{
  
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  if(dist %in% c("mbbefd", "MBBEFD"))
  {
    print("bllii")
  }else if(dist %in% c("oiunif", "oistpareto", "oibeta", "oigbeta")) #one inflated distr
  {
    p1 <- mean(x == 1)
    xneq1 <- x[x != 1]
    uplolist <- list(upper=Inf, lower=0)
    if(dist == "oistpareto")
    {
      start <- list(a=1)               
    }else if(dist == "oistpareto")
    {
      n <- length(xneq1)
      m <- mean(xneq1)
      v <- (n - 1)/n*var(xneq1)
      aux <- m*(1-m)/v - 1
      start <- list(shape1=m*aux, shape2=(1-m)*aux)
      
    }else if(dist == "oigbeta")
    {
      start <- list(shape0=1, shape1=1, shape2=1)
    }
      
    f1 <- fitdist(xneq1, distr=dist, method=method, start=start, fix.arg=list(p1=p1), 
                  lower=uplolist$lower, upper=uplolist$upper, ...)
    
  }else if(dist == "oiunif")
  {
    f1 <- fitdist(xneq1, distr=dist, method=method, start=list(p1=p1),
                  lower=0, upper=1, ...)
  }
}

