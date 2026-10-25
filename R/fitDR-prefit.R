# File overview: helpers of fitDR() (see fitDR.R):
# - prefitDR.mle(): starting values for the MBBEFD and GB1 fits;
# - d<dist>1/2/3(): densities on transformed (unconstrained) parameters, used
#   by the pre-fit; they must be exported because fitdistrplus::mledist()
#   looks up the density "d" + distr by name;
# - fitDR.selectbest(): choice of the best fit among the domains and Hessian;
# - tlmme.MBBEFD(): exact total-loss-moment-matching estimator of MBBEFD.

#search good starting values
#for MBBEFD, the parameter space is split in three domains (see fitDR()):
# D1: b > 1; D2: b < 1, gb < 1; D3: b < 1, gb > 1
# Pre-fit: maximize the likelihood without constraints, on transformed
# parameters u such that the original parameters are Trans.*(u) and always lie in
# the right domain (see util-CInftransform.R). The optimum, mapped back to the
# original parameters, is used as starting value of the constrained fit.
# The result is a list of starting values (one per domain), with NA when a
# pre-fit fails; fitDR() then keeps its default starting values.
prefitDR.mle <- function(x, dist, control=list(trace=0, REPORT=1, maxit=100), ...)
{
  if(dist == "mbbefd")
  {

    # MBBEFD(a,b): one pre-fit per domain, all starting from u = (0, 0)
    prefit1 <- mledist(x, distr="mbbefd1", optim.method="BFGS", start=list(a=0, b=0),
                       silent=TRUE, control=control, ...)
    prefit2 <- mledist(x, distr="mbbefd2", optim.method="BFGS", start=list(a=0, b=0),
                       silent=TRUE, control=control, ...)
    prefit3 <- mledist(x, distr="mbbefd3", optim.method="BFGS", start=list(a=0, b=0),
                       silent=TRUE, control=control, ...)

    # map the optimum back to the original parameters
    if(prefit1$convergence %in% 0:1) #either successful or reached the iter limit
    {
      initpar1 <- c(Trans.m10(prefit1$estimate["a"]), Trans.1Inf(prefit1$estimate["b"]))
    }else
    {
      initpar1 <- c(NA, NA)
    }

    if(prefit2$convergence %in% 0:1)
    {
      initpar2 <- c(Trans.0Inf(prefit2$estimate["a"]), Trans.01(prefit2$estimate["b"]))
    }else
    {
      initpar2 <- c(NA, NA)
    }

    if(prefit3$convergence %in% 0:1)
    {
      initpar3 <- c(Trans.mInfm1(prefit3$estimate["a"]), Trans.01(prefit3$estimate["b"]))
    }else
    {
      initpar3 <- c(NA, NA)
    }

    list(initpar1, initpar2, initpar3)

  }else if(dist == "MBBEFD")
  {
    #no constraint for the new param
    # MBBEFD(g,b), domain D1: g = 1 + exp(u1), b = 1 + exp(u2), no constraint needed
    prefit1 <- mledist(x, distr="MBBEFD1", optim.method="BFGS", start=list(g=0, b=0),
                       silent=TRUE, control=control, ...)

    # In D2 and D3, g = 1 + exp(u1) and b = logistic(u2) is not enough: gb < 1 (D2)
    # or gb > 1 (D3) is a linear constraint on (u1, u2), handled by stats::constrOptim().
    # constrOptim2 adapts constrOptim() to the custom.optim interface of mledist().
    constrOptim2 <- function(par, fn, gr=NULL, ui, ci, ...)
      constrOptim(theta=unlist(par), f=fn, grad=gr, ui=ui, ci=ci, ...)
    #constraint is  g+b < 0 <=> -g-b > 0 (new param)
    prefit2 <- mledist(x, distr="MBBEFD2", optim.method="BFGS", custom.optim=constrOptim2,
                       start=list(g=-1, b=0), ui = cbind(-1, -1), ci = 0,
                       silent=TRUE, control=control, ...)
    #constraint is gb > 1 <=> g+b > 0 (new param)
    prefit3 <- mledist(x, distr="MBBEFD3", optim.method="BFGS", custom.optim=constrOptim2,
                       start=list(g=1, b=0), ui = cbind(1, 1), ci = 0,
                       silent=TRUE, control=control, ...)
    # map the optimum back to the original parameters
    if(prefit1$convergence %in% 0:1) #either successful or reached the iter limit
    {
      initpar1 <- c(Trans.1Inf(prefit1$estimate["g"]), Trans.1Inf(prefit1$estimate["b"]))
    }else
    {
      initpar1 <- c(NA, NA)
    }

    if(prefit2$convergence %in% 0:1)
    {
      initpar2 <- c(Trans.1Inf(prefit2$estimate["g"]), Trans.01(prefit2$estimate["b"]))
    }else
    {
      initpar2 <- c(NA, NA)
    }

    if(prefit3$convergence %in% 0:1)
    {
      initpar3 <- c(Trans.1Inf(prefit3$estimate["g"]), Trans.01(prefit3$estimate["b"]))
    }else
    {
      initpar3 <- c(NA, NA)
    }

    list(initpar1, initpar2, initpar3)

  # GB1 (continuous part of oigbeta): the three shape parameters are positive,
  # shape = exp(u); total losses are removed first
  }else if(dist == "oigbeta")
  {
    x <- x[x != 1]
    prefit <- mledist(x, distr="gbeta1", optim.method="BFGS", silent=TRUE, control=control,
                  start=list(shape0=0, shape1=0, shape2=0))
    if(prefit$convergence %in% 0:1)
    {
      initpar <- Trans.0Inf(prefit$estimate)
    }else
    {
      initpar <- c(NA, NA)
    }
    initpar
  }else
    stop("wrong distribution in prefitDR.mle().")
}

# Densities on transformed parameters (u1, u2): each one maps (u1, u2) to the
# parameters of the corresponding domain and calls the original density.
#Transformed distribution (internal use)

#MBBEFD(a,b)
#domain : (a,b) in (-1, 0) x (1, +Inf)
dmbbefd1 <- function(x, a, b, log=FALSE)
  dmbbefd(x, Trans.m10(a), Trans.1Inf(b), log=log)
#domain : (a,b) in (0, +Inf) x (0, 1)
dmbbefd2 <- function(x, a, b, log=FALSE)
  dmbbefd(x, Trans.0Inf(a), Trans.01(b), log=log)
#domain : (a,b) in (-Inf, -1) x (0, 1)
dmbbefd3 <- function(x, a, b, log=FALSE)
  dmbbefd(x, Trans.mInfm1(a), Trans.01(b), log=log)

#MBBEFD(g,b)
#domain : (g,b) in (1, +Inf) x (1, +Inf)
#with gb > 1 (old param) : always verified in param
dMBBEFD1 <- function(x, g, b, log=FALSE)
  dMBBEFD(x, Trans.1Inf(g), Trans.1Inf(b), log=log)
#domain : (g,b) in (1, +Inf) x (0, 1)
# with gb < 1 (old param) : g < -b in new param
dMBBEFD2 <- function(x, g, b, log=FALSE)
  dMBBEFD(x, Trans.1Inf(g), Trans.01(b), log=log)
#domain : (g,b) in (1, +Inf) x (0, 1)
# with gb > 1 (old param) : g > -b in new param
dMBBEFD3 <- function(x, g, b, log=FALSE)
  dMBBEFD(x, Trans.1Inf(g), Trans.01(b), log=log)


#OI-GB1
dgbeta1 <- function(x, shape0, shape1, shape2, log=FALSE)
  dgbeta(x, Trans.0Inf(shape0), Trans.0Inf(shape1), Trans.0Inf(shape2), log=log)


# Keep the converged fit with the highest log-likelihood and compute the Hessian
# of -LL/n at the optimum (used by fitDR.addcomp() for the standard errors).
# When every domain fails, the first fit (with NA estimates) is returned.
#select the best MLE among several parameter domains and compute the Hessian of -LL/n
#fitlist: list of mledist() outputs (convergence == 100 means failure)
#hessian.fun: optional function returning the analytical Hessian of -LL/n at par
fitDR.selectbest <- function(fitlist, x, minusLL, grLL, hessian.fun=NULL)
{
  ok <- vapply(fitlist, function(f)
    !is.null(f$convergence) && f$convergence != 100 && isTRUE(is.finite(f$loglik)),
    logical(1))
  if(!any(ok))
    return(list(fit=fitlist[[1]], hessian=NULL))
  ll <- vapply(fitlist, function(f) if(isTRUE(is.finite(f$loglik))) f$loglik else -Inf, numeric(1))
  ll[!ok] <- -Inf
  f1 <- fitlist[[which.max(ll)]]
  if(!is.null(hessian.fun))
  {
    fHess <- try(hessian.fun(f1$estimate), silent=TRUE)
  }else
  {
    #steps relative to the parameter values: the default absolute step (1e-3) of optimHess() 
    #is too large for parameters close to a boundary or with small values
    ndeps <- pmax(abs(f1$estimate), 1e-6) * 1e-4
    fHess <- try(optimHess(par=f1$estimate, fn=minusLL, gr=grLL, obs=x, 
                           control=list(ndeps=ndeps)), silent=TRUE)
  }
  if(inherits(fHess, "try-error") || any(!is.finite(fHess)))
    fHess <- NULL
  list(fit=f1, hessian=fHess)
}

#exact total-loss-moment-matching estimator for MBBEFD(g,b), see Bernegger (1997), Section 4.1:
#g = 1/p where p = P(X = 1) is the empirical total loss, and
#b solves E[X] = mean(x), as E[X] is decreasing in b from 1 (b -> 0) to 1/g (b -> +Inf).
tlmme.MBBEFD <- function(x)
{
  # the system has a solution if and only if 0 < p < mean < 1 (the mean of an
  # MBBEFD distribution with P(X = 1) = p lies in (p, 1))
  p <- etl(x, na.rm=TRUE)
  m <- mean(x, na.rm=TRUE)
  if(!(p > 0 && p < 1 && m > p && m < 1))
  {
    warning("the total-loss-moment-matching estimator requires 0 < etl(x) < mean(x) < 1.")
    return(list(estimate=c(g=NA, b=NA), convergence=100))
  }
  g <- 1/p
  # Solve on log(b) so that the search covers b in (0, +Inf); extendInt widens
  # the initial interval as needed, using the fact that the function decreases.
  meandiff <- function(logb) mMBBEFD(1, g, exp(logb)) - m
  res <- try(uniroot(meandiff, interval=c(-5, 5), extendInt="downX",
                     tol=sqrt(.Machine$double.eps)*1e-4, maxiter=1000), silent=TRUE)
  if(inherits(res, "try-error") || !is.finite(res$root))
    return(list(estimate=c(g=NA, b=NA), convergence=100))
  list(estimate=c(g=g, b=exp(res$root)), convergence=0)
}
