bootDR <- function(f, bootmethod="param", niter=1001, silent=TRUE,
                   parallel = c("no", "snow", "multicore"), ncpus)
{
  if (!inherits(f, "DR"))
    stop("Use only with 'DR' objects")
  if (niter<10) 
    stop("niter must be an integer above 10")
  bootmethod <- match.arg(bootmethod, c("param", "nonparam"))
  
  parallel <- match.arg(parallel, c("no", "snow", "multicore"))
  if (parallel == "multicore" & .Platform$OS.type == "windows")
  {
    parallel <- "snow"
    warning("As the multicore option is not supported on Windows it was replaced by snow")
  }
  if ((parallel == "snow" | parallel == "multicore") & missing(ncpus)) 
    stop("You have to specify the number of available processors to parallelize 
             the bootstrap")
  
  #simulate bootstrap data
  if (bootmethod == "param") 
  { # parametric bootstrap
    rdistname <- paste("r", f$distname, sep="")
    if (!exists(rdistname, mode="function"))
      stop(paste("The ", rdistname, " function must be defined"))
    rdata <- do.call(rdistname, c(list(n=niter*f$n), as.list(f$estimate)))
    dim(rdata) <- c(f$n, niter)
  }else 
  { # non parametric bootstrap
    rdata <- sample(f$data, size=niter*f$n, replace=TRUE)
    dim(rdata) <- c(f$n, niter)
  }
  
  #compute bootstrap estimates
  start <- as.list(f$estimate) #a named vector is no longer accepted as starting values.
  
  func <- function(iter) 
  {
    res <- try(do.call(fitDR, list(x=rdata[, iter], dist=f$distname, start=start,
                                   optim.method=f$optim.method, control=f$control)
                       ), silent=silent)
    
    if(class(res)[1] == "try-error")
      return(c(rep(NA, length(start)), 100))
    else
      return(c(res$estimate, res$convergence))
  }
  
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  options(warn=ifelse(silent, -1, 0), show.error.messages=!silent)
  
  # parallel or sequential computation
  if (parallel != "no") 
  {
    if (parallel == "snow") type <- "PSOCK"
    else if (parallel == "multicore") type <- "FORK"
    clus <- parallel::makeCluster(ncpus, type = type)
    parallel::clusterEvalQ(clus, {library(mbbefd)})
    resboot <- parallel::parSapply(clus, 1:niter, func)
    parallel::stopCluster(clus)
  }
  else
  {
    resboot <- sapply(1:niter, func)
  }
  
  options(warn=owarn, show.error.messages=oerr)
  
  rownames(resboot) <- c(names(f$estimate), "convergence")
  if (length(resboot[, 1])>2) 
  {
    estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
    bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm=TRUE))
    colnames(bootCI) <- c("Median", "2.5%", "97.5%")
  }
  else {
    estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
    names(estim) <- names(f$estimate)
    bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm=TRUE)) 
    names(bootCI) <- c("Median", "2.5%", "97.5%") 
  } 
  
  # code of convergence of the optimization function for each iteration
  converg <- t(resboot)[, length(resboot[, 1])]
  
  res <- structure(list(estim=estim, converg=converg, 
                        method=bootmethod, nbboot=niter, CI=bootCI, fitpart=f), 
                   class=c("bootDR", "bootdist"))
  res    
}