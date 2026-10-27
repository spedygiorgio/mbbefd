# File overview: builds the object returned by fitDR(). The object mimics the
# output of fitdistrplus::fitdist() (same component names and order), so that
# the fitdist methods (summary, plot, cdfcomp, gofstat, bootdist, ...) can be
# used on fitDR objects.
# Arguments: x the data, theta the estimated parameters, hessian the Hessian of
# the average negative log-likelihood -LL/n at theta (NULL if not available),
# dist the distribution name, method "mle" or "tlmme", convergence the code of
# the optimizer.
fitDR.addcomp <- function(x, theta, hessian = NULL, dist, method, convergence = 0, vcov = NULL,
                          optim.method = "default", control = list())
{
  #components will be
  #"estimate", "method", "sd", "cor", "vcov", "loglik", "aic", "bic", "n", "data",
  #"distname", "fix.arg", "fix.arg.fun", "dots", "convergence", "discrete", "weights"
  #"optim.method", "control"

  f1 <- list(estimate = theta, weights = NULL, dots = NULL, fix.arg = NULL, fix.arg.fun = NULL)
  #other fitdist components
  f1$convergence <- convergence
  f1$method <- method
  f1$n <- length(x)
  f1$data <- x
  f1$distname <- dist
  f1$discrete <- FALSE
  f1$optim.method <- optim.method
  f1$control <- control
  npar <- length(theta)

  # log-likelihood and information criteria AIC = -2 LL + 2 k, BIC = -2 LL + log(n) k
  #gof statistics
  if (any(is.na(theta)))
  {
    f1$loglik <- f1$aic <- f1$bic <- NA
  } else
  {
    f1$loglik <- LLfunc(obs = x, theta = theta, dist = dist)
    f1$aic <- -2 * f1$loglik + 2 * npar
    f1$bic <- -2 * f1$loglik + log(f1$n) * npar
  }

  # covariance matrix of the estimates (only for maximum likelihood)
  #one-inflated uniform distribution
  if (dist %in% "oiunif")
  {
    stop("do not need to call fitDR.addcomp()")
  } else if (dist %in% c("oibeta", "oistpareto", "oigbeta"))
  {
    #one-inflated distribution with at least two parameters
    if (method == "mle")
    {
      # For one-inflated distributions, p1 is estimated by the proportion of total
      # losses, independently of the shape parameters, which are estimated on the
      # observations below 1 only. Hence the covariance matrix is block diagonal:
      # - shape parameters: inverse Hessian of -LL/n divided by the number of
      #   observations below 1;
      # - p1: binomial variance p1 (1 - p1) / n.
      p1 <- f1$estimate["p1"]
      if (is.null(hessian))
      {
        f1$vcov <- f1$sd <- f1$cor <- NA
      } else if (all(!is.na(hessian)) && qr(hessian)$rank == NCOL(hessian))
      {
        subvcov <- solve(hessian) / sum(x < 1)
        f1$vcov <- rbind(cbind(as.matrix(subvcov), rep(0, npar - 1)),
                         c(rep(0, npar - 1), p1 * (1 - p1) / f1$n))
        dimnames(f1$vcov) <- list(names(f1$estimate), names(f1$estimate))
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
      } else
      {
        f1$vcov <- f1$sd <- f1$cor <- NA
      }
    } else
    {
      f1$vcov <- f1$sd <- f1$cor <- NA
    }
  } else
  {
    # MBBEFD: the asymptotic covariance matrix is the inverse of the observed
    # information, i.e. solve(Hessian of -LL/n) / n.
    #non one-inflated distributions: mbbefd / MBBEFD
    if (method == "mle")
    {
      if (is.null(hessian))
        f1$vcov <- f1$sd <- f1$cor <- NA
      else if (all(!is.na(hessian)) && qr(hessian)$rank == NCOL(hessian))
      {
        f1$vcov <- solve(hessian)
        f1$vcov <- f1$vcov / f1$n
        dimnames(f1$vcov) <- list(names(f1$estimate), names(f1$estimate))
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
      } else
      {
        f1$vcov <- f1$sd <- f1$cor <- NA
      }
    } else
    {
      f1$vcov <- f1$sd <- f1$cor <- NA
    }

  }
  #output
  f1
}
