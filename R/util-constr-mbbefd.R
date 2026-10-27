# File overview: inequality constraints for the constrained maximum likelihood
# of the MBBEFD distribution in fitDR(). The optimizer alabama::constrOptim.nl()
# requires every constraint to be written as h(theta) > 0, where h returns a
# vector (argument hin) and its jacobian matrix d h / d theta (argument hin.jac,
# one row per constraint, one column per parameter).
#
# The MBBEFD parameter space is split into three domains, each fitted
# separately because the formulas change across the boundaries b = 1 and gb = 1:
#   D1: b > 1                (in (a,b): -1 < a < 0, b > 1)
#   D2: b < 1 and gb < 1     (in (a,b): a > 0, 0 < b < 1)
#   D3: b < 1 and gb > 1     (in (a,b): a < -1, 0 < b < 1)
# The signature (x, fix.arg, obs, ddistnam) matches the way fitdistrplus::mledist()
# calls the constraint functions: only x (the parameter vector) is used.

#to meet the standard 'fn' argument and specific name arguments,

#constraint function for MBBEFD(a,b)
# whole (a,b) domain (not used by fitDR)
constrmbbefd <- function(x, fix.arg, obs, ddistnam)
{
  res <- x[1] * (1 - x[2]) #a*(1-b) >= 0
  names(res) <- NULL
  res
}
#domain : (a,b) in (-1, 0) x (1, +Inf)
# domain D1 in (a,b)
constrmbbefd1 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1] + 1, -x[1], x[2] - 1, x[1] * (1 - x[2])) #-1 < a < 0, b > 1, a*(1-b) >= 0
  names(res) <- NULL
  res
}
# jacobian of constrmbbefd1: row i is the gradient of the i-th constraint
constrmbbefd1jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1, ] <- c(1, 0)      #d(a+1)
  j[2, ] <- c(-1, 0)     #d(-a)
  j[3, ] <- c(0, 1)      #d(b-1)
  j[4, ] <- c(1 - x[2], -x[1]) #d(a*(1-b))
  dimnames(j) <- NULL
  j
}

#domain : (a,b) in (0, +Inf) x (0, 1)
# domain D2 in (a,b)
constrmbbefd2 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1], x[2], 1 - x[2], x[1] * (1 - x[2])) #0 < a , 0 < b < 1, a*(1-b) >= 0
  names(res) <- NULL
  res
}
# jacobian of constrmbbefd2
constrmbbefd2jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1, ] <- c(1, 0)      #d(a)
  j[2, ] <- c(0, 1)      #d(b)
  j[3, ] <- c(0, -1)     #d(1-b)
  j[4, ] <- c(1 - x[2], -x[1]) #d(a*(1-b))
  dimnames(j) <- NULL
  j
}

#domain : (a,b) in (-Inf, -1) x (0, 1), i.e. b < 1 and gb > 1 in the 2nd param.
# domain D3 in (a,b)
constrmbbefd3 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(-1 - x[1], x[2], 1 - x[2]) #a < -1, 0 < b < 1
  names(res) <- NULL
  res
}
# jacobian of constrmbbefd3
constrmbbefd3jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 3, 2)
  j[1, ] <- c(-1, 0)
  j[2, ] <- c(0, 1)
  j[3, ] <- c(0, -1)
  dimnames(j) <- NULL
  j
}


#constraint function for MBBEFD(g,b)
# whole (g,b) domain (not used by fitDR)
constrMBBEFD <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1] - 1, x[2]) #g >= 1, b > 0
  names(res) <- NULL
  res
}
#domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
# domain D1 in (g,b)
constrMBBEFD1 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1] - 1, x[2] - 1, x[1] * x[2] - 1) #g > 1, b > 1, gb > 1
  names(res) <- NULL
  res
}
# jacobian of constrMBBEFD1
constrMBBEFD1jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 3, 2)
  j[1, ] <- c(1, 0)
  j[2, ] <- c(0, 1)
  j[3, ] <- c(x[2], x[1])
  dimnames(j) <- NULL
  j
}

#domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
# domain D2 in (g,b)
constrMBBEFD2 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1] - 1, 1 - x[2], x[2], 1 - x[1] * x[2]) #g > 1, 1 > b > 0, gb < 1
  names(res) <- NULL
  res
}
# jacobian of constrMBBEFD2
constrMBBEFD2jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1, ] <- c(1, 0)
  j[2, ] <- c(0, -1)
  j[3, ] <- c(0, 1)
  j[4, ] <- c(-x[2], -x[1])
  dimnames(j) <- NULL
  j
}

#domain : (g,b) in (1, +Inf) x (0, 1) with gb > 1 (decreasing density, see Gatti & Wuthrich, 2024)
# domain D3 in (g,b)
constrMBBEFD3 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1] - 1, 1 - x[2], x[1] * x[2] - 1) #g > 1, b < 1, gb > 1 (implies b > 0)
  names(res) <- NULL
  res
}
# jacobian of constrMBBEFD3
constrMBBEFD3jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 3, 2)
  j[1, ] <- c(1, 0)
  j[2, ] <- c(0, -1)
  j[3, ] <- c(x[2], x[1])
  dimnames(j) <- NULL
  j
}
