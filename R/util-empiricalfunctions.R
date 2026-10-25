# File overview: functions computed from a sample of destruction rates:
# the empirical exposure curve (class "eecf", with print/summary/plot/lines
# methods), the empirical total loss probability etl() and the Theil index.
######################################################################
# Empirical sample-based functions
######################################################################


#empirical exposure curve
#in the spirit of ecdf()
# Empirical exposure curve G_n(d) = mean(min(x_i, d)) / mean(x_i).
# The returned object is a function of d (like ecdf()); the sorted data and the
# values of G_n at the data points are stored in its environment for plotting.
eecf <- function(x)
{
  x <- sort(x) # drops NAs
  n <- length(x)
  if(n < 1) stop("'x' must have 1 or more non-missing values")
  
  # G_n at the sorted data points x_(k): (x_(1) + ... + x_(k) + (n-k) x_(k)) / n / mean(x)
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

# Print the call and the first/last data values.
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

# Summary statistics of the data behind the empirical exposure curve.
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

# Plot the empirical exposure curve as a continuous piecewise linear curve
# starting at (0,0), with the data points, the line y = 1 and the diagonal
# y = x (exposure curve of a total loss only distribution).
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

# add the empirical exposure curve to an existing plot
lines.eecf <- function(x, ...)
  plot(x, add=TRUE, ...)
  


#total loss
# Empirical total loss: proportion of observations equal to 1. It estimates
# P(X = 1), i.e. 1/g for MBBEFD(g,b) and p1 for one-inflated distributions.
etl <- function(x, na.rm = FALSE)
  mean(x == 1, na.rm = na.rm, trim=0)

###################
#internal functions

#Theil index, see package ineq for other income index (e.g. Gini coefficient)
# Empirical Theil index mean((x/m) log(x/m)), m = mean(x); used in fitDR() to
# get a starting value of the power parameter of the GB1 distribution.
Theil.emp <- function(x, na.rm = FALSE) 
  mean(x/mean(x, na.rm = na.rm, trim=0)*log(x/mean(x, na.rm = na.rm, trim=0)), 
       na.rm = na.rm, trim=0)
