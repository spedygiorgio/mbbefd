# File overview: bootstrap of destruction rate fits. bootDR() refits the model
# on many resampled datasets and returns the distribution of the estimates,
# in the same format as fitdistrplus::bootdist() (so that its methods work).
# Arguments: f a fitDR object, bootmethod "param" (simulate from the fitted
# distribution) or "nonparam" (resample the data with replacement), niter the
# number of bootstrap samples, silent whether to hide errors and warnings of the
# individual fits, parallel/ncpus to distribute the refits over several cores.
bootDR <- function(f, bootmethod = "param", niter = 1001, silent = TRUE,
                   parallel = c("no", "snow", "multicore"), ncpus)
{
  # check the arguments
  if (!inherits(f, "fitDR"))
    stop("Use only with 'fitDR' objects")
  if (niter < 10)
    stop("niter must be an integer above 10")
  bootmethod <- match.arg(bootmethod, c("param", "nonparam"))

  parallel <- match.arg(parallel, c("no", "snow", "multicore"))
  if (parallel == "multicore" && .Platform$OS.type == "windows")
  {
    parallel <- "snow"
    warning("As the multicore option is not supported on Windows it was replaced by snow")
  }
  if ((parallel == "snow" || parallel == "multicore") && missing(ncpus))
    stop("You have to specify the number of available processors to parallelize 
             the bootstrap")

  # Generate all the bootstrap samples at once, as the columns of a
  # matrix with f$n rows and niter columns.
  #simulate bootstrap data
  if (bootmethod == "param")
  { # parametric bootstrap
    rdistname <- paste("r", f$distname, sep = "")
    if (!exists(rdistname, mode = "function"))
      stop(paste("The ", rdistname, " function must be defined"))
    rdata <- do.call(rdistname, c(list(n = niter * f$n), as.list(f$estimate)))
    dim(rdata) <- c(f$n, niter)
  } else
  { # non parametric bootstrap
    rdata <- sample(f$data, size = niter * f$n, replace = TRUE)
    dim(rdata) <- c(f$n, niter)
  }

  #compute bootstrap estimates
  start <- as.list(f$estimate) #a named vector is no longer accepted as starting values.

  # refit the model on the iter-th bootstrap sample; return the estimates and the
  # convergence code (100 and NA estimates when the fit fails)
  func <- function(iter)
  {
    res <- try(do.call(fitDR, list(x = rdata[, iter], dist = f$distname, start = start,
                                   optim.method = f$optim.method, control = f$control)
                       ), silent = silent)

    if (inherits(res, "try-error"))
      c(rep(NA, length(start)), 100)
    else
      c(res$estimate, res$convergence)
  }

  # temporarily silence warnings and error messages of the individual fits
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  options(warn = ifelse(silent, -1, 0), show.error.messages = !silent)
  #restore the user options even if the bootstrap loop fails (e.g. an interruption)
  on.exit(options(warn = owarn, show.error.messages = oerr), add = TRUE)

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

  options(warn = owarn, show.error.messages = oerr)

  # resboot has one column per bootstrap sample: the parameters, then the
  # convergence code. Summarize it by the median and the 2.5% / 97.5% quantiles
  # (percentile bootstrap confidence intervals).
  rownames(resboot) <- c(names(f$estimate), "convergence")
  if (length(resboot[, 1]) > 2)
  {
    estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
    bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm = TRUE),
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm = TRUE),
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm = TRUE))
    colnames(bootCI) <- c("Median", "2.5%", "97.5%")
  }
  # one-parameter case (resboot has a single parameter row)
  else {
    estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
    names(estim) <- names(f$estimate)
    bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm = TRUE),
                quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm = TRUE),
                quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm = TRUE))
    names(bootCI) <- c("Median", "2.5%", "97.5%")
  }

  # code of convergence of the optimization function for each iteration
  converg <- t(resboot)[, length(resboot[, 1])]

  res <- structure(list(estim = estim, converg = converg,
                        method = bootmethod, nbboot = niter, CI = bootCI, fitpart = f),
                   class = c("bootDR", "bootdist"))
  res
}



# Covariance matrix of the bootstrap estimates.
vcov.bootDR <- function(object, ...)
{
  stopifnot(inherits(object, "bootDR"))
  cov(object$estim, ...)
}
