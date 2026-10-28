# File overview: user-level d, p, q, r, ec (exposure curve), m (raw moment) and
# tl (total loss) functions of the one-inflated uniform distribution: with
# probability p1 the loss is a total loss (X = 1), otherwise X is uniform on [0,1).
# Every function is a thin wrapper around the generic one-inflated functions of
# distr-1infl-distrib.R, applied to the base R uniform functions with min=0, max=1.
#d, p, q, r function for one-inflated uniform distribution


# density (point mass p1 at x = 1)
doiunif <- function(x, p1, log = FALSE)
{
  doifun(x = x, dfun = dunif, p1 = p1, log = log, min = 0, max = 1)
}

# distribution function
poiunif <- function(q, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q = q, pfun = punif, p1 = p1, lower.tail = lower.tail, log.p = log.p, min = 0, max = 1)
}

# quantile function
qoiunif <- function(p, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p = p, qfun = qunif, p1 = p1, lower.tail = lower.tail, log.p = log.p, min = 0, max = 1)
}

# random generation
roiunif <- function(n, p1)
{
  roifun(n = n, rfun = runif, p1 = p1, min = 0, max = 1)
}


# exposure curve
ecoiunif <- function(x, p1)
{
  ecoifun(x = x, ecfun = ecunif, mfun = munif, p1 = p1, min = 0, max = 1)
}

# raw moment of a given order
moiunif <- function(order, p1)
{
  moifun(order = order, mfun = munif, p1 = p1, min = 0, max = 1)
}

# total loss probability P(X = 1), i.e. p1
tloiunif <- function(p1)
{
  tloifun(p1 = p1)
}
