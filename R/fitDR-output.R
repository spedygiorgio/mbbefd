fitDR.addcomp <- function(x, theta, hessian=NULL, dist, method, convergence=0, vcov=NULL)
{
  #components will be
  #"estimate", "method", "sd", "cor", "vcov", "loglik", "aic", "bic", "n", "data", 
  #"distname", "fix.arg", "fix.arg.fun", "dots", "convergence", "discrete", "weights"
  
  
  f1 <- list(estimate=theta, weights=NULL, dots=NULL, fix.arg=NULL, fix.arg.fun=NULL, dots=NULL)
  #other fitdist components
  f1$convergence <- convergence
  f1$method <- method
  f1$n <- length(x)
  f1$data <- x
  f1$distname <- dist
  f1$discrete <- FALSE
  
  
  #gof stat
  f1$loglik <- LLfunc(obs=x, theta=theta, dist=dist)
  npar <- length(theta)
  f1$aic <- -2*f1$loglik+2*npar
  f1$bic <- -2*f1$loglik+log(f1$n)*npar
  
  #non one-inflated distributions
  if(!dist %in% c("oiunif", "oibeta", "oistpareto", "oigbeta"))
  {
    #computes Hessian of -log Lik at estimate values
    if(all(!is.na(theta)) && method == "mle" && !is.null(hessian))
    {
      if(all(!is.na(hessian)) && qr(hessian)$rank == NCOL(hessian)){
        f1$vcov <- solve(hessian)
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
      }else
        f1$vcov <- f1$sd <- f1$cor <- NA
    }else{
      f1$vcov <- f1$sd <- f1$cor <- NA
    }
  }else #one-inflated distributions
  {
    p1 <- f1$estimate["p1"]
    if(method == "mle")
    {  
      if(all(!is.na(vcov)))
      {
        f1$vcov <- rbind(cbind(as.matrix(vcov), rep(0, npar-1)),
                         c(rep(0, npar-1), p1*(1-p1)))
        f1$sd <- sqrt(diag(f1$vcov))
        f1$cor <- cov2cor(f1$vcov)
      }
    }else
    {
      f1$vcov <- rbind(cbind(matrix(NA, npar-1, npar-1), rep(0, npar-1)),
                       c(rep(0, npar-1), p1*(1-p1)))
      f1$sd <- f1$cor <- NA
    }
    dimnames(f1$vcov) <- list(names(f1$estimate), names(f1$estimate))
  }
  #output
  f1
}