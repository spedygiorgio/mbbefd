
fitDR <- function(x, dist, method="mle", form.arg=NULL, start=NULL, ...)
{
  if(any(x < 0 | x > 1))
    stop("Values outside [0,1] are not supported in fitDR.")
  method <- match.arg(method, c("mle", "mme", "qme", "mge"))
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  
  if(dist %in% c("mbbefd", "MBBEFD"))
  {
    stop("not yet implemented")
  }else if(dist == "oiunif")
  {
    if(is.null(start))
      start=list(p1=etl(x))
    
    #print(llfunc(x, start$p1, dist))
    f1 <- fitdist(x, distr=dist, method=method, start=start,
                  lower=0, upper=1, ..., optim.method="Brent") #, control=list(trace=6, REPORT=1)
    class(f1) <- c("DR", class(f1))
    
  }else if(dist %in% c("oistpareto", "oibeta", "oigbeta")) #one-inflated distr
  {
    p1 <- etl(x)
    xneq1 <- x[x != 1]
    distneq1 <- substr(dist, 3, nchar(dist))
    
    #print(dist)
    #print(distneq1)
    
    uplolist <- list(upper=Inf, lower=0)
    if(is.null(start))
    {
      if(distneq1 == "stpareto")
      {
        start <- list(a=1)               
      }else if(distneq1 == "beta")
      {
        n <- length(xneq1)
        m <- mean(xneq1)
        v <- (n - 1)/n*var(xneq1)
        aux <- m*(1-m)/v - 1
        start <- list(shape1=m*aux, shape2=(1-m)*aux)
        
      }else if(distneq1 == "gbeta")
      {
        start <- list(shape0=1, shape1=1, shape2=1)
      }else
        stop("wrong non-inflated distribution.")
    }else
    {
      if(distneq1 == "stpareto")
      {
        start <- start["a"]               
      }else if(distneq1 == "beta")
      {
        start <- start[c("shape1", "shape2")]
      }else if(distneq1 == "gbeta")
      {
        start <- start[c("shape0", "shape1", "shape2")]
      }else
        stop("wrong non-inflated distribution.")
    }
    #print(start)
      
    if(method == "mle")
    {
      f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
                  lower=uplolist$lower, upper=uplolist$upper, ...)
      if(f1$convergence != 0)
      {
         stop("error in convergence when fitting data.")
      }else
      {
        f1$estimate <- c(f1$estimate, p1=p1) 
        f1$n <- length(x)
        f1$distname <- dist
        f1$data <- x
        
        #gof stat
        f1$loglik <- llfunc(obs=x, theta=f1$estimate, dist=dist)
        npar <- length(f1$estimate)
        f1$aic <- -2*f1$loglik+2*npar
        f1$bic <- -2*f1$loglik+log(f1$n)*npar
        
        f1$vcov <- rbind(cbind(as.matrix(f1$vcov), rep(0, npar-1)), 
                         c(rep(0, npar-1), p1*(1-p1)))
        dimnames(f1$vcov) <- list(names(f1$estimate), names(f1$estimate))
        
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
        class(f1) <- c("DR", class(f1))
      } 
      
    }else
    {
      stop("not yet implemented.")
    }
    
  }else
    stop("Unknown distribution for destruction rate models.")
  
  f1
}


llfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  ddist <- paste0("d", dist)
  sum(log(do.call(ddist, c(list(obs), as.list(theta)) ) ) )
}


