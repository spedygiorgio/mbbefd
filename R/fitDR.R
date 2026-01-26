
fitDR <- function(x, dist, method="mle", start=NULL, optim.method="default", control=list(trace=0), ...)
{
  
  if(any(is.na(x)))
    x <- x[!is.na(x)]
  if(length(x) <= 0)
    stop("There is no data in 'x' argument (possibly after removing NA).")
  if(any(x < 0 | x > 1))
    stop("Values outside [0,1] are not supported in fitDR.")
  method <- match.arg(method, c("mle", "tlmme"))
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  
  #make control parameters as in optim()
  con4optim <- list(trace = 0, fnscale = 1, maxit = 100L, abstol = -Inf, 
                    reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
                    gamma = 2, REPORT = 1, warn.1d.NelderMead = TRUE, type = 1, 
                    lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  idx <- names(control) %in% names(con4optim)
  con4optim[names(control)[idx]] <- control[names(control)[idx]]
  con4optim$trace <- max(con4optim$trace -1, 0) #decrease trace level
  
  #make control parameters as in constrOptim.nl()
  con4constrOptim.nl <- list(mu0 = 0.01, sig0 = 10, eps = 1e-07, itmax = 50, 
                             method = "BFGS", trace = FALSE, NMinit = FALSE)
  idx <- names(control) %in% names(con4constrOptim.nl)
  con4constrOptim.nl[names(control)[idx]] <- control[names(control)[idx]]
  con4constrOptim.nl$trace <- max(con4constrOptim.nl$trace -1, 0)
  
  if(is.null(control$trace))
    control$trace <- 0
  if(control$trace > 0)
    cat("Distribution is", dist, "\n")
  
  if(dist == "mbbefd")
  {
    initparmbbefd <- list(list(a=Trans.m10(0), b=Trans.1Inf(0)), 
                          list(a=Trans.0Inf(0), b=Trans.01(0)))
    if(control$trace > 0)
    {
      cat("\tinitial parameter value try\n")
      print(sapply(initparmbbefd, unlist))
    }
    prefit <- prefitDR.mle(x, "mbbefd")
    #domain : (a,b) in (-1, 0) x (1, +Inf)
    if(all(!is.na(prefit[[1]])))
    {
      initparmbbefd[[1]] <- prefit[[1]]
      if(initparmbbefd[[1]]["b"] == 1)
        initparmbbefd[[1]]["b"] <- 2
    }
    #domain : (a,b) in (0, +Inf) x (0, 1)
    if(all(!is.na(prefit[[2]])))
    {
      initparmbbefd[[2]] <- prefit[[2]]
      if(initparmbbefd[[2]]["b"] == 1)
        initparmbbefd[[2]]["b"] <- 1/2
    }
    if(!is.list(initparmbbefd[[1]]))
      initparmbbefd[[1]] <- as.list(initparmbbefd[[1]])
    if(!is.list(initparmbbefd[[2]]))
      initparmbbefd[[2]] <- as.list(initparmbbefd[[2]])
    if(control$trace > 0)
    {
      cat("\tinitial parameter value after prefit()\n")
      print(sapply(initparmbbefd, unlist))
    }
    
    if(method == "mle")
    {
      
      #wrap -LL to match the call by fitdist
      minusLL <- function(x, fix.arg, obs, ddistnam) 
        -LLfunc(obs=obs, theta=x, dist="mbbefd")/length(obs)
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) 
        -grLLfunc(obs=obs, theta=x, dist="mbbefd")/length(obs)
      
      #domain : (a,b) in (-1, 0) x (1, +Inf)
      alabama1 <- mledist(x, distr="mbbefd", start=initparmbbefd[[1]], custom.optim= constrOptim.nl, 
                          hin=constrmbbefd1, control.outer=con4constrOptim.nl, gradient=grLL, 
                          calcvcov=FALSE, silent = con4constrOptim.nl$trace == 0, ...)
      #domain : (a,b) in (0, +Inf) x (0, 1)
      alabama2 <- mledist(x, distr="mbbefd", start=initparmbbefd[[2]], custom.optim= constrOptim.nl, 
                          hin=constrmbbefd2, control.outer=con4constrOptim.nl, gradient=grLL, 
                          calcvcov=FALSE, silent = con4constrOptim.nl$trace == 0, ...)
      fHess <- NULL
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
        
        fHess <- try(optimHess(par=f1$estimate, fn=minusLL, obs=x, gr=grLL), silent=TRUE)
        if(inherits(fHess, "try-error"))
          fHess <- NULL
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=fHess, vcov=NULL,
                          dist="mbbefd", method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mmbbefd(1, par[1], par[2]) - mean(obs))^2 + (tlmbbefd(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparmbbefd[[1]]), fn=DIFF2, hin=constrmbbefd1, 
                                 hin.jac=constrmbbefd1jac, obs=x, control.outer=con4constrOptim.nl)
      alabama2 <- constrOptim.nl(unlist(initparmbbefd[[2]]), fn=DIFF2, hin=constrmbbefd2, 
                                 hin.jac=constrmbbefd2jac, obs=x, control.outer=con4constrOptim.nl)
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
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, vcov=NULL,
                          dist="mbbefd", method="tlmme", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after optim()\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")
  }else if(dist == "MBBEFD")
  {
    #starting values
    g <- 1/etl(x, na.rm=TRUE)
    if(is.infinite(g))
      g <- 2
    initparMBBEFD <- list(list(g=g, b=Trans.1Inf(0)), 
                          list(g=g, b=1/(2*g)))
    if(control$trace > 0)
    {
      cat("\tinitial parameter value first try\n")
      print(sapply(initparMBBEFD, unlist))
    }
    #try to improve initial value for b
    phalf <- mean(x <= 1/2)
    b <- Re(polyroot(c(phalf, (1-g)*(1-phalf), - 1 +g*(1-phalf))))
    if(any(b > 0 | b < 1))
      initparMBBEFD[[2]]["b"] <- max(b[b > 0 | b < 1])
    else if(any(b > 1))
      initparMBBEFD[[2]]["b"] <- max(b[b > 1])
    if(control$trace > 0)
    {
      cat("\tinitial parameter value second try\n")
      print(sapply(initparMBBEFD, unlist))
    }
    prefit <- prefitDR.mle(x, "MBBEFD")
    if(all(!is.na(prefit[[1]])))
    { 
      initparMBBEFD[[1]] <- prefit[[1]]
      if(initparMBBEFD[[1]]["b"] == 1)
        initparMBBEFD[[1]]["b"] <- 2
    }
    if(all(!is.na(prefit[[2]])))
    {  
      initparMBBEFD[[2]] <- prefit[[2]]
      if(initparMBBEFD[[2]]["b"] == 1)
        initparMBBEFD[[2]]["b"] <- 1/2
    }
    if(!is.list(initparMBBEFD[[1]]))
      initparMBBEFD[[1]] <- as.list(initparMBBEFD[[1]])
    if(!is.list(initparMBBEFD[[2]]))
      initparMBBEFD[[2]] <- as.list(initparMBBEFD[[2]])
    if(control$trace > 0)
    {
      cat("\tinitial parameter value after prefit()\n")
      print(sapply(initparMBBEFD, unlist))
    }
    if(method == "mle")
    {
      #wrap -LL to match the call by fitdist
      minusLL <- function(x, fix.arg, obs, ddistnam) 
        -LLfunc(obs=obs, theta=x, dist="MBBEFD")/length(obs)
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) 
        -grLLfunc(obs=obs, theta=x, dist="MBBEFD")/length(obs)
      
      #domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
      alabama1 <- mledist(x, distr="MBBEFD", start=initparMBBEFD[[1]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD1, 
                          control.outer=con4constrOptim.nl, hin.jac=constrMBBEFD1jac, 
                          silent = con4constrOptim.nl$trace == 0, calcvcov=FALSE, gradient=grLL)
      #domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
      alabama2 <- mledist(x, distr="MBBEFD", start=initparMBBEFD[[2]], 
                          custom.optim= constrOptim.nl, hin=constrMBBEFD2, 
                          control.outer=con4constrOptim.nl, hin.jac=constrMBBEFD2jac, 
                          silent = con4constrOptim.nl$trace == 0, calcvcov=FALSE, gradient=grLL)
      
      fHess <- NULL
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
        fHess <- try(optimHess(par=f1$estimate, fn=minusLL, obs=x, gr=grLL), silent=FALSE)
        if(inherits(fHess, "try-error"))
          fHess <- NULL
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=fHess, vcov=NULL,
                          dist="MBBEFD", method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }
    }else if(method == "tlmme")
    {
      DIFF2 <- function(par, obs) 
      {
        (mMBBEFD(1, par[1], par[2]) - mean(obs))^2 + (tlMBBEFD(par[1], par[2]) - etl(obs))^2
      }
      alabama1 <- constrOptim.nl(unlist(initparMBBEFD[[1]]), fn=DIFF2, hin=constrMBBEFD1, 
                                 obs=x, control.outer=con4constrOptim.nl)
      alabama2 <- constrOptim.nl(unlist(initparMBBEFD[[2]]), fn=DIFF2, hin=constrMBBEFD2, 
                                 obs=x, control.outer=con4constrOptim.nl)
      
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
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, vcov=NULL,
                          dist="MBBEFD", method="tlmme", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after optim()\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")  
    
  }else if(dist == "oiunif")
  {
    if(is.null(start))
      start <- list(p1=etl(x))
    if(control$trace > 0)
    {
      cat("\tinitial parameter value try\n")
      print(unlist(start))
    }
    if(method %in% c("mle", "tlmme"))
    {
      if(method == "tlmme")
        method <- "mle"
      if(optim.method == "default")
        optim.method <- "Brent"
      f1 <- fitdist(x, distr=dist, method=method, start=start, calcvcov=TRUE, 
                    lower=0, upper=1, control=con4optim, ..., optim.method=optim.method) 
      f1$optim.method <- optim.method
      f1$control <- NULL
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after fitdist()\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")  
    
  }else if(dist %in% c("oistpareto", "oibeta", "oigbeta")) #one-inflated distr
  {
    p1 <- etl(x, na.rm=TRUE)
    xneq1 <- x[x != 1]
    distneq1 <- substr(dist, 3, nchar(dist))
    
    if(is.null(start))
    {
      if(distneq1 == "stpareto")
      {
        start <- list(a=1)               
      }else if(distneq1 == "beta")
      {
        n <- length(xneq1)
        if(n > 1)
        {
          m <- mean(xneq1, na.rm=TRUE)
          v <- (n - 1)/n*var(xneq1, na.rm=TRUE)
          aux <- m*(1-m)/v - 1
          start <- list(shape1=m*aux, shape2=(1-m)*aux)
        }else
          start <- list(shape1=1, shape2=1) #i.e. uniform distribution
        
      }else if(distneq1 == "gbeta")
      {
        n <- length(xneq1)
        if(n > 1)
        {
          shape00 <- optimize(function(z) (Theil.emp(x, na.rm=TRUE) - Theil.theo.shape0(z, obs=x))^2, 
                              lower=0.01, upper=100)$minimum
          otherpar00 <- mmedist(x^shape00, "beta", method="mme", calcvcov = FALSE, control=con4optim)
          start <- c(list(shape0=shape00), as.list(otherpar00$estimate))
        }else
          start <- list(shape0=1, shape1=1, shape2=1) #i.e. uniform distribution
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
    if(control$trace > 0)
    {
      cat("\tinitial parameter value try\n")
      print(unlist(start))
    }  
    if(method == "mle")
    {
      uplolist <- list(upper=Inf, lower=0)
      #check the initial value
      loglik0 <- LLfunc(xneq1, unlist(start), distneq1)
      
      if(is.infinite(loglik0))
        warning("initial value of the log-likelihood is infinite.")
      
      #improve initial parameters for GB1
      if(distneq1 == "gbeta")
      {
        if(length(xneq1) > 1)
        {
          prefit <- prefitDR.mle(xneq1, "oigbeta")
          
          if(all(!is.na(prefit)))
            start <- as.list(prefit)
          if(control$trace > 0)
          {
            cat("\tinitial parameter value after prefit()\n")
            print(unlist(start))
          }
        }
        if(optim.method == "default")
          optim.method <- "BFGS"
        f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, control=con4optim,
                      optim.method=optim.method, calcvcov=TRUE, ...)
        
      }else
      {
        if(optim.method == "default")
          optim.method <- "L-BFGS-B"
        f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
                      lower=uplolist$lower, upper=uplolist$upper, control=con4optim,
                      optim.method=optim.method, calcvcov=TRUE, ...)
      }
      f1$estimate <- c(f1$estimate, "p1"=p1) 
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, vcov=NULL,
                          dist=dist, method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after fitdist()\n")
        print(f1$estimate)
      }  
      
    }else if(method == "tlmme")
    {
      start <- c(start, list(p1=p1))
      if(control$trace > 0)
      {
        cat("\tinitial parameter value try\n")
        print(unlist(start))
      }
      npar <- length(start)
      #param should be positive, and p1 <= 1
      uplolist <- list(upper=c(rep(Inf, npar-1), 1), lower=0)
      
      DIFF2 <- function(par, obs) 
      {
        #true 
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
      res <- optim(par=unlist(start), fn=DIFF2, obs=x, method="L-BFGS-B", control=con4optim,
                   lower=uplolist$lower, upper=uplolist$upper, ...)
      if(res$convergence > 0 && res$value > sqrt(.Machine$double.eps))
        f1 <- list(estimate=rep(NA, npar), convergence=100, hessian=NULL)
      else
      {
        f1 <- list(estimate=res$par, convergence=0, hessian=NULL)
      }
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, vcov=NULL,
                          dist=dist, method="tlmme", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after optim()\n")
        print(f1$estimate)
      }
    }else
    {
      stop("not yet implemented.")
    }
    
  }else
    stop("Unknown distribution for destruction rate models.")
  
  #reorder components as a fitdist object
  f1 <- f1[c("estimate", "method", "sd", "cor", "vcov", "loglik", "aic", "bic", "n", "data", 
             "distname", "fix.arg", "fix.arg.fun", "dots", "convergence", "discrete", "weights",
             "optim.method", "control")]
  class(f1) <- c("DR", "fitdist")
  f1
}

