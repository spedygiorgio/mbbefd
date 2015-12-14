
fitDR <- function(x, dist, method="mle", start=NULL, ...)
{
  if(any(is.na(x)))
    x <- x[!is.na(x)]
  if(any(x < 0 | x > 1))
    stop("Values outside [0,1] are not supported in fitDR.")
  method <- match.arg(method, c("mle", "tlmme"))
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  
  cat("blii\n")
  if(dist == "mbbefd")
  {
    initparmbbefd <- list(list(a=-1/2, b=2), list(a=2, b=1/2))
    
    if(method == "mle")
    {
      
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) -grLLfunc(obs=obs, theta=x, dist="mbbefd")
      
      #domain : (a,b) in (-1, 0) x (1, +Inf)
      alabama1 <- fitdist(x, distr="mbbefd", start=initparmbbefd[[1]], 
                        custom.optim= constrOptim.nl, hin=constrmbbefd1, method="mle",
                        control.outer=list(trace= FALSE), gr=grLL)
      #domain : (a,b) in (0, +Inf) x (0, 1)
      alabama2 <- fitdist(x, distr="mbbefd", start=initparmbbefd[[2]], 
                        custom.optim= constrOptim.nl, hin=constrmbbefd2, method="mle",
                        control.outer=list(trace= FALSE), gr=grLL)
      if(alabama1$convergence == 100 && alabama2$convergence == 100)
        f1 <- alabama1
      else if(alabama1$convergence == 100 && alabama2$convergence != 100) 
        f1 <- alabama2
      else if(alabama1$convergence != 100 && alabama2$convergence == 100) 
        f1 <- alabama1
      else
      {
        if(alabama1$loglik > alabama2$loglik)
          f1 <- alabama1
        else
          f1 <- alabama2  
        #computes Hessian of -LL at estimate values
        f1hess <- -heLLfunc(obs=x, theta=f1$estimate, dist="mbbefd")
        
        if(all(!is.na(f1hess)) && qr(f1hess)$rank == NCOL(f1hess)){
          f1$vcov <- solve(f1hess)
          f1$sd <- sqrt(diag(f1$vcov))
          f1$cor <- cov2cor(f1$vcov)
        }#otherwise it is already at NA from fitdist
      }
      
      class(f1) <- c("DR", class(f1))
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mmbbefd(1, par[1], par[2]) - mean(obs))^2 + (tlmbbefd(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparmbbefd[[1]]), fn=DIFF2, hin=constrmbbefd1, obs=x, control.outer=list(trace=FALSE))
      alabama2 <- constrOptim.nl(unlist(initparmbbefd[[2]]), fn=DIFF2, hin=constrmbbefd2, obs=x, control.outer=list(trace=FALSE))
      
      if(alabama1$convergence > 0 && alabama2$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else if(alabama1$convergence > 0 && alabama2$convergence == 0) 
        f1 <- list(estimate=alabama2$par, convergence=0)
      else if(alabama1$convergence == 0 && alabama2$convergence > 0) 
        f1 <- list(estimate=alabama1$par, convergence=0)
      else
      {
        if(alabama1$value < alabama2$value)
          f1 <- list(estimate=alabama1$par, convergence=0)
        else
          f1 <- list(estimate=alabama2$par, convergence=0)
      }
      f1$method <- "tlmme"
      f1$data <- x
      f1$n <- length(x)
      f1$distname <- "mbbefd"
      f1$fix.arg <- f1$fix.arg.fun <- f1$dots <- f1$weights <- NULL
      f1$discrete <- FALSE
      f1$aic <- f1$bic <- f1$loglik <- NA
      f1$sd <- f1$vcov <- f1$cor <- NA
      class(f1) <- c("DR", "fitdist")
      
    }else
      stop("not yet implemented")
  }else if(dist == "MBBEFD")
  {
    xneq1 <- x[x != 1]
    g <- 1/etl(x, na.rm=TRUE)
    
    phalf <- mean(x <= 1/2)
    b <- Re(polyroot(c(phalf, (1-g)*(1-phalf), - 1 +g*(1-phalf))))
    if(any(b > 0 | b < 1))
      b <- max(b[b > 0 & b < 1])
    else
      b <- 1/2
    
    initparMBBEFD <- list(list(g=1/etl(x), b=2), list(g=1/etl(x), b=etl(x)/2))
    
    #cat("g", g, "b", b, "\n")
    if(method == "mle")
    {
      
      #domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
      alabama1 <- fitdist(x, distr="MBBEFD", start=initparMBBEFD[[1]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD1, method="mle",
                          control.outer=list(trace= FALSE), hin.jac=constrMBBEFD1jac, silent=FALSE)
      #domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
      alabama2 <- fitdist(x, distr="MBBEFD", start=initparMBBEFD[[2]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD2, method="mle",
                          control.outer=list(trace= FALSE), hin.jac=constrMBBEFD2jac, silent=FALSE)
      
      print(summary(alabama1))
      print(summary(alabama2))
      if(alabama1$convergence == 100 && alabama2$convergence == 100)
        f1 <- alabama1
      else if(alabama1$convergence == 100 && alabama2$convergence != 100) 
        f1 <- alabama2
      else if(alabama1$convergence != 100 && alabama2$convergence == 100) 
        f1 <- alabama1
      else
      {
        if(alabama1$loglik > alabama2$loglik)
          f1 <- alabama1
        else
          f1 <- alabama2 
        #gof stat
        f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
        npar <- length(f1$estimate)
        f1$aic <- -2*f1$loglik+2*npar
        f1$bic <- -2*f1$loglik+log(f1$n)*npar
      }
      
      class(f1) <- c("DR", class(f1))
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mMBBEFD(1, par[1], par[2]) - mean(obs))^2 + (tlMBBEFD(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparMBBEFD[[1]]), fn=DIFF2, hin=constrMBBEFD1, obs=x, control.outer=list(trace=FALSE))
      alabama2 <- constrOptim.nl(unlist(initparMBBEFD[[2]]), fn=DIFF2, hin=constrMBBEFD2, obs=x, control.outer=list(trace=FALSE))
      
      if(alabama1$convergence > 0 && alabama2$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else if(alabama1$convergence > 0 && alabama2$convergence == 0) 
        f1 <- list(estimate=alabama2$par, convergence=0)
      else if(alabama1$convergence == 0 && alabama2$convergence > 0) 
        f1 <- list(estimate=alabama1$par, convergence=0)
      else
      {
        if(alabama1$value < alabama2$value)
          f1 <- list(estimate=alabama1$par, convergence=0)
        else
          f1 <- list(estimate=alabama2$par, convergence=0)
      }
      f1$method <- "tlmme"
      f1$data <- x
      f1$n <- length(x)
      f1$distname <- "MBBEFD"
      f1$fix.arg <- f1$fix.arg.fun <- f1$dots <- f1$weights <- NULL
      f1$discrete <- FALSE
      #gof stat
      f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
      npar <- length(f1$estimate)
      f1$aic <- -2*f1$loglik+2*npar
      f1$bic <- -2*f1$loglik+log(f1$n)*npar
      
      f1$sd <- f1$vcov <- f1$cor <- NA
      class(f1) <- c("DR", "fitdist")
      
    }else
      stop("not yet implemented")  
      
  }else if(dist == "oiunif")
  {
    if(is.null(start))
      start <- list(p1=etl(x))
    
    #print(LLfunc(x, start$p1, dist))
    
    if(method %in% c("mle", "tlmme"))
    {
      if(method == "tlmme")
        method <- "mle"
      f1 <- fitdist(x, distr=dist, method=method, start=start,
                  lower=0, upper=1, ..., optim.method="Brent") #, control=list(trace=6, REPORT=1)
    
      class(f1) <- c("DR", class(f1))
    }else
      stop("not yet implemented")  
    
  }else if(dist %in% c("oistpareto", "oibeta", "oigbeta")) #one-inflated distr
  {
    cat("zog\n")
    
    p1 <- etl(x, na.rm=TRUE)
    xneq1 <- x[x != 1]
    distneq1 <- substr(dist, 3, nchar(dist))
    
    #print(dist)
    #print(distneq1)
    cat("tau\n")
    
    uplolist <- list(upper=Inf, lower=0)
    if(is.null(start))
    {
      if(distneq1 == "stpareto")
      {
        start <- list(a=1)               
      }else if(distneq1 == "beta")
      {
        n <- length(xneq1)
        m <- mean(xneq1, na.rm=TRUE)
        v <- (n - 1)/n*var(xneq1, na.rm=TRUE)
        aux <- m*(1-m)/v - 1
        start <- list(shape1=m*aux, shape2=(1-m)*aux)
        
      }else if(distneq1 == "gbeta")
      {
        shape00 <- optimize(function(z) (Theil.emp(x, na.rm=TRUE) - Theil.theo.shape0(z, obs=x))^2, lower=0.01, upper=20)$minimum
        start <- c(list(shape0=shape00), as.list(fitdist(x^shape00, "beta", method="mme")$estimate))
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
      cat(distneq1, "\n")
      print(start)
      f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
                  lower=uplolist$lower, upper=uplolist$upper, ...)
      cat("blii\n")
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
        f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
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
      print(f1)
      
    }else if(method == "tlmme")
    {
      start <- c(start, list(p1=p1))
      npar <- length(start)
      
      DIFF2 <- function(par, obs) 
      {
        PX1 <- do.call(paste0("tl", dist), as.list(par))
        EX <- do.call(paste0("m", dist), as.list(c(order=1, par)))
        if(npar <= 2)
          return( (EX - mean(obs))^2 + (PX1 - etl(obs))^2 )
        
        if(npar >= 3)
          EX2 <- do.call(paste0("m", dist), as.list(c(order=2, par)))
        if(npar >= 4)
          EX3 <- do.call(paste0("m", dist), as.list(c(order=3, par)))
        
        if(npar == 3)
          return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (PX1 - etl(obs))^2 )
        else if(npar == 4)
          return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (EX3 - mean(obs^3))^2 + (PX1 - etl(obs))^2 )
        else
          stop("not implemented")
      }
          
      
      res <- optim(par=unlist(start), fn=DIFF2, obs=x, method="L-BFGS-B", lower=uplolist$lower, upper=uplolist$upper)
      
      if(res$convergence > 0)
        f1 <- list(estimate=NA, convergence=100)
      else
      {
        f1 <- list(estimate=res$par, convergence=0)
      }
      f1$method <- "tlmme"
      f1$data <- x
      f1$n <- length(x)
      f1$distname <- dist
      f1$fix.arg <- f1$fix.arg.fun <- f1$dots <- f1$weights <- NULL
      f1$discrete <- FALSE
      #gof stat
      f1$loglik <- LLfunc(obs=x, theta=f1$estimate, dist=dist)
      
      f1$aic <- -2*f1$loglik+2*npar
      f1$bic <- -2*f1$loglik+log(f1$n)*npar
      
      f1$sd <- f1$vcov <- f1$cor <- NA
      class(f1) <- c("DR", "fitdist")
      
    }else
    {
      stop("not yet implemented.")
    }
    
  }else
    stop("Unknown distribution for destruction rate models.")
  
  f1
}

#likelihood function
LLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  ddist <- paste0("d", dist)
  sum(log(do.call(ddist, c(list(obs), as.list(theta)) ) ) )
}

#gradient of the likelihood function
grLLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("mbbefd", "MBBEFD")) 
  if(dist == "mbbefd")
  {
    g1 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, (b-1)/(a+1)/(a+b), (2*a+1)/(a*(a+1)) - 2/(a+b^x))
    }
    g2 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, a/(b*(a+b)), -x/b+1/(b*log(b))+2*a*x/(b*(a+b^x)))
    }
    c(sum(sapply(obs, g1, theta=theta)), sum(sapply(obs, g2, theta=theta)))
  }else
  {
    stop("not yet implemented.")
  }
}

#Hessian of the likelihood function
heLLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("mbbefd", "MBBEFD")) 
  if(dist == "mbbefd")
  {
    h11 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2-1/(a+1)^2, 2/(a+b^x)^2-1/a^2-1/(a+1)^2)
    }
    h21 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2, 2*x*b^(x-1)/(a+b^x)^2)
    }
    h22 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2-1/b^2, 
             x/b^2-(log(b)+1)/(b^2*log(b)^2)-2*a*x/(b^2*(a+b^x))-2*a*x^2*b^x/(b^2*(a+b^x)^2))
    }
    rbind(c(sum(sapply(obs, h11, theta=theta)), sum(sapply(obs, h21, theta=theta))),
    c(sum(sapply(obs, h21, theta=theta)), sum(sapply(obs, h22, theta=theta))))
  }else
  {
    stop("not yet implemented.")
  }
}

#to meet the standard 'fn' argument and specific name arguments, 

#constraint function for MBBEFD(a,b)
constrmbbefd <- function(x, fix.arg, obs, ddistnam)
{
  x[1]*(1-x[2]) #a*(1-b) >= 0
}
constrmbbefd1 <- function(x, fix.arg, obs, ddistnam)
{
  c(x[1]+1, -x[1], x[2]-1, x[1]*(1-x[2])) #-1 < a < 0, b > 1, a*(1-b) >= 0
}
constrmbbefd2 <- function(x, fix.arg, obs, ddistnam)
{
  c(x[1], x[1], 1-x[2], x[1]*(1-x[2])) #0 < a , 0 < b < 1, a*(1-b) >= 0
}

#constraint function for MBBEFD(g,b)
constrMBBEFD <- function(x, fix.arg, obs, ddistnam)
{
  c(x[1]-1, x[2]) #g >= 1, b > 0
}
constrMBBEFD1 <- function(x, fix.arg, obs, ddistnam)
{
  c(x[1]-1, x[2]-1, x[1]*x[2]-1) #g > 1, b > 1, gb > 1
}
constrMBBEFD1jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 3, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, 1)
  j[3,] <- c(x[2], x[1])
  j
}

constrMBBEFD2 <- function(x, fix.arg, obs, ddistnam)
{
  c(x[1]-1, 1-x[2], x[2], 1-x[1]*x[2]) #g > 1, 1 > b > 0, gb < 1
}
constrMBBEFD2jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, -1)
  j[3,] <- c(0, 1)
  j[4,] <- c(-x[2], -x[1])
  j
}
constrMBBEFDb <- function(x, fix.arg, obs, ddistnam)
{
  x[1] #b > 0
}
