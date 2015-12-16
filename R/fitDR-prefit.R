
#search good starting values
prefitDR.mle <- function(x, dist, control=list(trace=0, REPORT=1, maxit=100), ...)
{
  if(dist == "mbbefd")
  {
    
    prefit1 <- mledist(x, "mbbefd1", optim.method="BFGS", start=list(a=0, b=0), 
                       silent=TRUE, control=control, ...)
    prefit2 <- mledist(x, "mbbefd2", optim.method="BFGS", start=list(a=0, b=0), 
                       silent=TRUE, control=control, ...)
    
    if(prefit1$convergence == 0)
    {
      initpar1 <- c(Trans.m10(prefit1$estimate["a"]), Trans.1Inf(prefit1$estimate["b"]))
    }else
    {
      initpar1 <- c(NA, NA)
    }
    
    if(prefit2$convergence == 0)
    {
      initpar2 <- c(Trans.0Inf(prefit2$estimate["a"]), Trans.01(prefit2$estimate["b"]))  
    }else
    {
      initpar2 <- c(NA, NA)  
    }
    
    list(initpar1, initpar2)
    
  }else if(dist == "MBBEFD")
  {
    prefit1 <- mledist(x, "MBBEFD1", optim.method="BFGS", start=list(g=0, b=0), 
                       silent=TRUE, control=control, ...)
    prefit2 <- mledist(x, "MBBEFD2", optim.method="BFGS", start=list(g=0, b=0), 
                       silent=TRUE, control=control, ...)
    print(str(prefit1))
    print(str(prefit2))
    if(prefit1$convergence == 0)
    {
      initpar1 <- c(Trans.1Inf(prefit1$estimate["g"]), Trans.1Inf(prefit1$estimate["b"]))
    }else
    {
      initpar1 <- c(NA, NA)
    }
    
    if(prefit2$convergence == 0)
    {
      initpar2 <- c(Trans.1Inf(prefit2$estimate["g"]), Trans.01(prefit2$estimate["b"]))  
    }else
    {
      initpar2 <- c(NA, NA)  
    }
    
    list(initpar1, initpar2)
    
  }else if(dist == "oigbeta")
  {
    x <- x[x != 1] 
    prefit <- mledist(x, dist="gbeta1", optim.method="BFGS", silent=TRUE, control=control, 
                  start=list(shape0=0, shape1=0, shape2=0)) 
    if(prefit$convergence == 0)
    {
      initpar <- Trans.0Inf(prefit2$estimate)
    }else
    {
      initpar <- c(NA, NA)  
    }
    initpar
  }
}

#Transformed distribution (internal use)

#MBBEFD(a,b)
#domain : (a,b) in (-1, 0) x (1, +Inf)
dmbbefd1 <- function(x, a, b, log=FALSE)
  dmbbefd(x, Trans.m10(a), Trans.1Inf(b), log=log)
#domain : (a,b) in (0, +Inf) x (0, 1)
dmbbefd2 <- function(x, a, b, log=FALSE)
  dmbbefd(x, Trans.0Inf(a), Trans.01(b), log=log)

#MBBEFD(a,b)
#domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
dMBBEFD1 <- function(x, g, b, log=FALSE)
  dMBBEFD(x, Trans.1Inf(g), Trans.1Inf(b), log=log)
#domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
dMBBEFD2 <- function(x, g, b, log=FALSE)
  dMBBEFD(x, Trans.1Inf(g), Trans.01(b), log=log)


#OI-GB1
dgbeta1 <- function(x, shape0, shape1, shape2, log=FALSE)
  dgbeta(x, Trans.0Inf(shape0), Trans.0Inf(shape1), Trans.0Inf(shape2), log=log)
