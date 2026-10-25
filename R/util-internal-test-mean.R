# File overview: three independent ways of computing the expectation E[X] of a
# one-inflated distribution on [0,1]. They are not exported: they are only used
# in the test scripts to check that the moment functions (m<dist>) are consistent
# with the density, the distribution function and the exposure curve.
#internal test function of expectation for one-inflated distribution

# E[X] from the density: integral of x f(x) on [0,1) plus 1 * P(X = 1),
# where P(X = 1) is the value returned by the "density" at x = 1 (point mass).
tmean1 <- function(doifun, ...)
  integrate(function(x) x*doifun(x, ...), 0, 1)$value + 1*doifun(x=1, ...)
# E[X] from the distribution function: for X in [0,1], E[X] = integral of 1 - F(x).
tmean2 <- function(poifun, ...)
  integrate(function(x) 1-poifun(x, ...), 0, 1)$value
# E[X] from the exposure curve G: G'(0) = 1/E[X], approximated by eps / G(eps).
tmean3 <- function(ecoifun, eps=1e-9, ...)
  eps/ecoifun(eps, ...)

