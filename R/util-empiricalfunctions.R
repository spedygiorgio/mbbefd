######################################################################
# Empirical sample-based functions
######################################################################


#empirical exposure curve
#in the spirit of ecdf()
eecf <- function(x)
{
  x <- sort(x) # drops NAs
  n <- length(x)
  if(n < 1) stop("'x' must have 1 or more non-missing values")
  
  Gx <- cumsum(x)/n + x*(n-(1:n))/n #numerator
  Gx <- Gx/mean(x) #denominator
  
  f <- function(d)
    mean(pmin(x, d))/mean(x)
  rval <- Vectorize(f, "d")
  
  class(rval) <- c("eecf", class(rval))
  assign("Gx", Gx, envir=environment(rval))
  assign("x", x, envir=environment(rval))
  assign("nobs", n, envir=environment(rval)) # e.g. to reconstruct rank(x)
  attr(rval, "call") <- sys.call()
  rval
}

print.eecf <- function (x, digits = getOption("digits") - 2L, ...)
{
  numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")
  cat("Empirical Exposure Curve Function \nCall: ")
  print(attr(x, "call"), ...)
  xx <- environment(x)$"x"
  n <- environment(x)$"nobs"
  
  i1 <- 1L:min(3L,n)
  i2 <- if(n >= 4L) max(4L, n-1L):n else integer()
  cat(" x[1:",n,"] = ", numform(xx[i1]),
      if(n>3L) ", ", if(n>5L) " ..., ", numform(xx[i2]), "\n", sep = "")
  invisible(x)
}

summary.eecf <- function(object, ...)
{
  xx <- environment(object)$"x"
  n <- environment(object)$"nobs"
  
  header <- paste("Empirical Exposure Curve Function:   ", n,
                  "unique values with summary\n")
  structure(summary(xx, digits = max(3, getOption("digits")-3), ...),
            header = header, class = "summary.eecf")
}

print.summary.eecf <- function(x, ...)
{
  cat(attr(x, "header"))
  y <- unclass(x); attr(y, "header") <- NULL
  print(y, ...)
  invisible(x)
}

plot.eecf <- function(x, ..., ylab="Gn(x)", do.points=TRUE, 
                      col.01line = "gray70", pch = 19, main=NULL,
                      ylim=NULL, add=FALSE)
{
  n <- environment(x)$"nobs"
  xx <- environment(x)$"x"
  yy <- environment(x)$"Gx"
  #add left-hand point
  if(length(xx) != n)
    stop("wrong x attribute")
  xx <- c(0, xx)
  yy <- c(0, yy)
  
  #from plot.stepfun called for plot.ecdf
  if(missing(main))
  main <- {
    cl <- attr(x,"call")
    deparse(if(!is.null(cl))cl else sys.call())
  }
  if(missing(ylim))
    ylim <- c(0, 1)
  
  #unlike ecdf, eecf is a continuous function
  if(!add)
  {  
    plot(xx, yy, type = "l", ylab = ylab, main=main, ylim=ylim, ...)
    if(do.points) points(xx[-1], yy[-1], pch = pch, ...)
    abline(h = 1, col = col.01line, lty = 2)
    abline(a = 0, b = 1, col = col.01line, lty = 2)
  }else
  {
    lines(xx, yy, type="l", ...)
    if(do.points) points(xx[-1], yy[-1], pch = pch, ...)
  }
  #terminates with invisible()
}

lines.eecf <- function(x, ...)
  plot(x, add=TRUE, ...)
  


#total loss
etl <- function(x, na.rm = FALSE)
  mean(x == 1, na.rm = na.rm, trim=0)

###################
#internal functions

#Theil index, see package ineq for other income index (e.g. Gini coefficient)
Theil.emp <- function(x, na.rm = FALSE) 
  mean(x/mean(x, na.rm = na.rm, trim=0)*log(x/mean(x, na.rm = na.rm, trim=0)), 
       na.rm = na.rm, trim=0)
