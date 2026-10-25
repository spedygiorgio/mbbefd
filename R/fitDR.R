# File overview: fitDR(), the main fitting function of the package.
# It fits a destruction rate model (a distribution on [0,1] with a point mass
# at 1 for total losses) to a sample x of destruction rates (loss / sum insured).
#
# Supported distributions (argument dist):
#   "mbbefd" / "MBBEFD"  MBBEFD distribution in its (a,b) / (g,b) parametrization
#   "oiunif", "oistpareto", "oibeta", "oigbeta"  one-inflated distributions
# Supported methods:
#   "mle"    maximum likelihood estimation;
#   "tlmme"  total-loss-moment matching: match the empirical total loss
#            probability etl(x) and the empirical moments.
#
# How the MBBEFD fit works (method "mle"):
# 1. the parameter space is split in three domains D1 (b > 1), D2 (b < 1, gb < 1)
#    and D3 (b < 1, gb > 1), because the parameter constraints differ;
# 2. starting values in each domain come from an unconstrained pre-fit on
#    transformed parameters (prefitDR.mle() in fitDR-prefit.R);
# 3. in each domain, -LL is minimized under inequality constraints with
#    alabama::constrOptim.nl(), called through fitdistrplus::mledist() and using
#    the analytical gradient (util-loglikfunc.R) and constraints
#    (util-constr-mbbefd.R);
# 4. each domain fit is refined by an unconstrained maximization in
#    (log(g-1), log(b)), where the log-likelihood is continuous across the domain
#    boundaries (the barrier of constrOptim.nl() stops short of them);
# 5. the fit with the highest log-likelihood is retained and its Hessian gives
#    the standard errors (fitDR.selectbest() in fitDR-prefit.R).
# For the one-inflated distributions, p1 is the proportion of total losses and
# the other parameters are fitted on the observations below 1.
# The result is built by fitDR.addcomp() (fitDR-output.R) as a "fitDR" object,
# which inherits from fitdistrplus "fitdist" objects.

fitDR <- function(x, dist, method="mle", start=NULL, optim.method="default", control=list(trace=0), ...)
{
  
  # Check the data: remove missing values, require at least one value in [0,1].
  if(any(is.na(x)))
    x <- x[!is.na(x)]
  if(length(x) <= 0)
    stop("There is no data in 'x' argument (possibly after removing NA).")
  if(any(x < 0 | x > 1))
    stop("Values outside [0,1] are not supported in fitDR.")
  method <- match.arg(method, c("mle", "tlmme"))
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD"))
  
  # The user control list is split into the options understood by optim() (used
  # by the pre-fits and the one-inflated fits) and those understood by
  # constrOptim.nl() (used by the MBBEFD fits). The trace level is decreased by
  # one for the optimizers, so that trace = 1 only prints fitDR's own messages.
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
  
  # Without total losses, the MLE of P(X = 1) is 0, i.e. on the boundary of the
  # parameter space (g -> +Inf): the fit is still computed but the user is warned.
  if(dist %in% c("mbbefd", "MBBEFD") && method == "mle" && !any(x == 1))
    warning("there is no total loss (x == 1) in the data: the estimated total loss probability goes to the boundary 0.")
  
  # ---------------------------------------------------------------------------
  # MBBEFD(a,b)
  # ---------------------------------------------------------------------------
  if(dist == "mbbefd")
  {
    #the parameter space is split in three domains
    #D1 : (a,b) in (-1, 0) x (1, +Inf)
    #D2 : (a,b) in (0, +Inf) x (0, 1)
    #D3 : (a,b) in (-Inf, -1) x (0, 1), i.e. b < 1 and gb > 1 in the 2nd parametrization
    # default starting values, one per domain (overwritten by the pre-fit below)
    initparmbbefd <- list(list(a=Trans.m10(0), b=Trans.1Inf(0)), 
                          list(a=Trans.0Inf(0), b=Trans.01(0)),
                          list(a=Trans.mInfm1(0), b=Trans.01(0)))
    if(control$trace > 0)
    {
      cat("\tinitial parameter value try\n")
      print(sapply(initparmbbefd, unlist))
    }
    if(method == "mle")
    {
      # pre-fit on transformed parameters to get better starting values
      prefit <- prefitDR.mle(x, "mbbefd", control=con4optim)
      #domain D1
      if(all(!is.na(prefit[[1]])))
      {
        initparmbbefd[[1]] <- prefit[[1]]
        if(initparmbbefd[[1]]["b"] == 1)
          initparmbbefd[[1]]["b"] <- 2
      }
      #domain D2 and D3
      for(i in 2:3)
        if(all(!is.na(prefit[[i]])))
        {
          initparmbbefd[[i]] <- prefit[[i]]
          if(initparmbbefd[[i]]["b"] == 1)
            initparmbbefd[[i]]["b"] <- 1/2
        }
      # mledist() requires the starting values as a named list
      for(i in 1:3)
        if(!is.list(initparmbbefd[[i]]))
          initparmbbefd[[i]] <- as.list(initparmbbefd[[i]])
      if(control$trace > 0)
      {
        cat("\tinitial parameter value after prefit()\n")
        print(sapply(initparmbbefd, unlist))
      }
      
      # -LL/n and its gradient, with the argument names expected by mledist()
      #wrap -LL to match the call by fitdist
      minusLL <- function(x, fix.arg, obs, ddistnam) 
        -LLfunc(obs=obs, theta=x, dist="mbbefd")/length(obs)
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) 
        -grLLfunc(obs=obs, theta=x, dist="mbbefd")/length(obs)
      
      #a failure in one domain is expected: do not print the error unless tracing
      oerr <- options(show.error.messages = con4constrOptim.nl$trace > 0)
      on.exit(options(oerr), add=TRUE)
      # one constrained fit per domain
      hinlist <- list(constrmbbefd1, constrmbbefd2, constrmbbefd3)
      hinjaclist <- list(constrmbbefd1jac, constrmbbefd2jac, constrmbbefd3jac)
      fitlist <- lapply(1:3, function(i)
        mledist(x, distr="mbbefd", start=initparmbbefd[[i]], custom.optim= constrOptim.nl, 
                hin=hinlist[[i]], hin.jac=hinjaclist[[i]], control.outer=con4constrOptim.nl, 
                gradient=grLL, calcvcov=FALSE, silent = con4constrOptim.nl$trace == 0, ...))
      options(oerr)
      # keep the best of the three fits and compute its Hessian, then build the output
      #analytical Hessian of -LL/n
      heLL <- function(par) -heLLfunc(obs=x, theta=par, dist="mbbefd")/length(x)
      # each domain fit is refined without constraints in (g,b), see fitDR.polish.mbbefd()
      best <- fitDR.selectbest(fitlist, x=x, minusLL=minusLL, grLL=grLL, hessian.fun=heLL,
                               polish=function(par) fitDR.polish.mbbefd(x, par))
      f1 <- best$fit
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=best$hessian, vcov=NULL,
                          dist="mbbefd", method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }
    }else if(method == "tlmme")
    {
      #exact solution of E[X] = mean(x) and P(X=1) = etl(x), see Bernegger (1997)
      # the exact solution is computed in (g,b) and converted to a = (g-1) b / (1 - gb)
      f1 <- tlmme.MBBEFD(x)
      gb <- f1$estimate
      if(f1$convergence == 0)
        f1$estimate <- c(a=unname((gb["g"]-1)*gb["b"]/(1-gb["g"]*gb["b"])), b=unname(gb["b"]))
      else
        f1$estimate <- c(a=NA, b=NA)
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=NULL, vcov=NULL,
                          dist="mbbefd", method="tlmme", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after tlmme\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")
  # ---------------------------------------------------------------------------
  # MBBEFD(g,b): same procedure as for MBBEFD(a,b)
  # ---------------------------------------------------------------------------
  }else if(dist == "MBBEFD")
  {
    #the parameter space is split in three domains
    #D1 : (g,b) in (1, +Inf) x (1, +Inf)
    #D2 : (g,b) in (1, +Inf) x (0, 1) with gb < 1
    #D3 : (g,b) in (1, +Inf) x (0, 1) with gb > 1
    if(method == "mle")
    {
      #starting values
      # g = 1 / P(X = 1) is estimated by the inverse of the empirical total loss
      g <- 1/etl(x, na.rm=TRUE)
      if(is.infinite(g))
        g <- 2
      # default starting values, one per domain (overwritten by the pre-fit below)
      initparMBBEFD <- list(list(g=g, b=Trans.1Inf(0)), 
                            list(g=g, b=1/(2*g)),
                            list(g=g, b=(1+1/g)/2))
      if(control$trace > 0)
      {
        cat("\tinitial parameter value first try\n")
        print(sapply(initparMBBEFD, unlist))
      }
      #try to improve initial value for b using F(1/2) = P(X <= 1/2), 
      #which is a polynomial equation of degree 2 in y = sqrt(b)
      phalf <- mean(x <= 1/2)
      y <- polyroot(c(phalf, (1-g)*(1-phalf), - 1 +g*(1-phalf)))
      y <- Re(y[abs(Im(y)) < 1e-8])
      b <- y[y > 0]^2
      if(any(b < 1/g))
        initparMBBEFD[[2]]["b"] <- max(b[b < 1/g])
      else if(any(b > 1))
        initparMBBEFD[[1]]["b"] <- max(b[b > 1])
      if(control$trace > 0)
      {
        cat("\tinitial parameter value second try\n")
        print(sapply(initparMBBEFD, unlist))
      }
      # pre-fit on transformed parameters to get better starting values
      prefit <- prefitDR.mle(x, "MBBEFD", control=con4optim)
      if(all(!is.na(prefit[[1]])))
      { 
        initparMBBEFD[[1]] <- prefit[[1]]
        if(initparMBBEFD[[1]]["b"] == 1)
          initparMBBEFD[[1]]["b"] <- 2
      }
      for(i in 2:3)
        if(all(!is.na(prefit[[i]])))
        {  
          initparMBBEFD[[i]] <- prefit[[i]]
          if(initparMBBEFD[[i]]["b"] == 1)
            initparMBBEFD[[i]]["b"] <- 1/2
        }
      # mledist() requires the starting values as a named list
      for(i in 1:3)
        if(!is.list(initparMBBEFD[[i]]))
          initparMBBEFD[[i]] <- as.list(initparMBBEFD[[i]])
      if(control$trace > 0)
      {
        cat("\tinitial parameter value after prefit()\n")
        print(sapply(initparMBBEFD, unlist))
      }
      
      # -LL/n and its gradient, with the argument names expected by mledist()
      #wrap -LL to match the call by fitdist
      minusLL <- function(x, fix.arg, obs, ddistnam) 
        -LLfunc(obs=obs, theta=x, dist="MBBEFD")/length(obs)
      #wrap gradient -LL to match the call by fitdist
      grLL <- function(x, fix.arg, obs, ddistnam) 
        -grLLfunc(obs=obs, theta=x, dist="MBBEFD")/length(obs)
      
      #a failure in one domain is expected: do not print the error unless tracing
      oerr <- options(show.error.messages = con4constrOptim.nl$trace > 0)
      on.exit(options(oerr), add=TRUE)
      # one constrained fit per domain
      hinlist <- list(constrMBBEFD1, constrMBBEFD2, constrMBBEFD3)
      hinjaclist <- list(constrMBBEFD1jac, constrMBBEFD2jac, constrMBBEFD3jac)
      fitlist <- lapply(1:3, function(i)
        mledist(x, distr="MBBEFD", start=initparMBBEFD[[i]], 
                custom.optim= constrOptim.nl, hin=hinlist[[i]], hin.jac=hinjaclist[[i]],
                control.outer=con4constrOptim.nl, gradient=grLL,
                silent = con4constrOptim.nl$trace == 0, calcvcov=FALSE, ...))
      options(oerr)
      # keep the best of the three fits, compute its Hessian (finite differences of
      # the analytical gradient), then build the output
      # each domain fit is refined without constraints, see fitDR.polish.MBBEFD()
      best <- fitDR.selectbest(fitlist, x=x, minusLL=minusLL, grLL=grLL,
                               polish=function(par) fitDR.polish.MBBEFD(x, par))
      f1 <- best$fit
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=best$hessian, vcov=NULL,
                          dist="MBBEFD", method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }
    }else if(method == "tlmme")
    {
      #exact solution of E[X] = mean(x) and P(X=1) = etl(x), see Bernegger (1997)
      f1 <- tlmme.MBBEFD(x)
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=NULL, vcov=NULL,
                          dist="MBBEFD", method="tlmme", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after tlmme\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")
    
  # ---------------------------------------------------------------------------
  # One-inflated uniform: the only parameter is p1; for both methods, the MLE
  # is the proportion of total losses, computed by fitdist() with Brent's method.
  # ---------------------------------------------------------------------------
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
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }
    }else
      stop("not yet implemented")  
    
  # ---------------------------------------------------------------------------
  # Other one-inflated distributions: p1 is estimated by etl(x), the parameters
  # of the continuous part (distneq1 = "stpareto", "beta" or "gbeta") are fitted
  # on the observations strictly below 1 (xneq1).
  # ---------------------------------------------------------------------------
  }else if(dist %in% c("oistpareto", "oibeta", "oigbeta")) #one-inflated distr
  {
    p1 <- etl(x, na.rm=TRUE)
    xneq1 <- x[x < 1]
    distneq1 <- substr(dist, 3, nchar(dist))
    
    # default starting values of the continuous part: a = 1 for stpareto, moment
    # estimates for beta, Theil index + moment estimates for GB1
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
    # user starting values: keep only the parameters of the continuous part
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
    # maximum likelihood on xneq1 with box constraints (parameters > 0)
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
          prefit <- prefitDR.mle(xneq1, "oigbeta", control=con4optim)
          
          loglik1 <- LLfunc(xneq1, prefit, distneq1)
          
          if(all(!is.na(prefit)) && loglik1 > loglik0)
            start <- as.list(prefit)
          if(control$trace > 0)
          {
            cat("\tinitial parameter value after prefit()\n")
            print(unlist(start))
          }
        }
        if(optim.method == "default")
          optim.method <- "L-BFGS-B"
        f1 <- mledist(xneq1, distr=distneq1, lower=uplolist$lower, upper=uplolist$upper, 
                      start=start, control=con4optim, optim.method=optim.method, calcvcov=TRUE, ...)
        #f1 <- fitdist(xneq1, distr=distneq1, method="mle", lower=uplolist$lower, upper=uplolist$upper, 
        #              start=start, control=con4optim,
        #              optim.method=optim.method, calcvcov=TRUE, ...)
        
      }else
      {
        if(optim.method == "default")
          optim.method <- "L-BFGS-B"
        f1 <- mledist(xneq1, distr=distneq1, lower=uplolist$lower, upper=uplolist$upper, 
                      start=start, control=con4optim, optim.method=optim.method, calcvcov=TRUE, ...)
        #f1 <- fitdist(xneq1, distr=distneq1, method="mle", start=start, 
        #              lower=uplolist$lower, upper=uplolist$upper, control=con4optim,
        #              optim.method=optim.method, calcvcov=TRUE, ...)
      }
      f1$estimate <- c(f1$estimate, "p1"=p1) 
      f1 <- fitDR.addcomp(x=x, theta=f1$estimate, hessian=f1$hessian, vcov=NULL,
                          dist=dist, method="mle", convergence=f1$convergence,
                          optim.method=optim.method, control=control)
      if(control$trace > 0)
      {
        cat("\toptimal parameter value after mledist()\n")
        print(f1$estimate)
      }  
      
    # match the total loss probability and the first raw moments (as many
    # equations as parameters) by least squares with optim(); fails (NA) when the
    # equations cannot be solved
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
  class(f1) <- c("fitDR", "fitdist")
  f1
}

