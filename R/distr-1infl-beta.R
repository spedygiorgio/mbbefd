# File overview: user-level d, p, q, r, ec (exposure curve), m (raw moment) and
# tl (total loss) functions of the one-inflated beta distribution: with
# probability p1 the loss is a total loss (X = 1), otherwise X is beta distributed
# on [0,1). Every function is a thin wrapper around the generic one-inflated
# functions of distr-1infl-distrib.R, applied to the base R beta functions.
# The non-centrality parameter ncp is only supported by d, p, q, r.
#d, p, q, r function for one-inflated beta distribution


# density (point mass p1 at x = 1)
doibeta <- function(x, shape1, shape2, p1, ncp = 0, log = FALSE)
{
  doifun(x = x, dfun = dbeta, p1 = p1, log = log, shape1 = shape1, shape2 = shape2, ncp = ncp)
}

# distribution function
poibeta <- function(q, shape1, shape2, p1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q = q, pfun = pbeta, p1 = p1, lower.tail = lower.tail, log.p = log.p,
         shape1 = shape1, shape2 = shape2, ncp = ncp)
}

# quantile function
qoibeta <- function(p, shape1, shape2, p1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p = p, qfun = qbeta, p1 = p1, lower.tail = lower.tail, log.p = log.p,
         shape1 = shape1, shape2 = shape2, ncp = ncp)
}

# random generation
roibeta <- function(n, shape1, shape2, p1, ncp = 0)
{
  roifun(n = n, rfun = rbeta, p1 = p1, shape1 = shape1, shape2 = shape2, ncp = ncp)
}

# exposure curve (central beta only)
ecoibeta <- function(x, shape1, shape2, p1, ncp = 0)
{
  if (ncp != 0)
    stop("not yet implemented.")
  ecoifun(x = x, ecfun = ecbeta, mfun = mbeta, p1 = p1, shape1 = shape1, shape2 = shape2)
}

# raw moment of a given order (central beta only)
moibeta <- function(order, shape1, shape2, p1, ncp = 0)
{
  if (ncp != 0)
    stop("not yet implemented.")
  moifun(order = order, mfun = mbeta, p1 = p1, shape1 = shape1, shape2 = shape2)
}


# total loss probability P(X = 1), i.e. p1
tloibeta <- function(shape1, shape2, p1, ncp = 0)
{
  tloifun(p1 = p1, shape1 = shape1, shape2 = shape2, ncp = ncp)
}
