# File overview: user-level d, p, q, r, ec (exposure curve), m (raw moment) and
# tl (total loss) functions of the one-inflated generalized beta of the first kind
# (GB1): with probability p1 the loss is a total loss (X = 1), otherwise X follows
# the GB1 distribution of distr-GB1.R (shape0 is the power parameter, shape1 and
# shape2 the beta parameters). Every function is a thin wrapper around the generic
# one-inflated functions of distr-1infl-distrib.R.
# (The comment below says "beta" but this file is for the generalized beta.)
#d, p, q, r function for one-inflated beta distribution


# density (point mass p1 at x = 1)
doigbeta <- function(x, shape0, shape1, shape2, p1, log = FALSE)
{
  doifun(x = x, dfun = dgbeta, p1 = p1, log = log, shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# distribution function
poigbeta <- function(q, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q = q, pfun = pgbeta, p1 = p1, lower.tail = lower.tail, log.p = log.p,
         shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# quantile function
qoigbeta <- function(p, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p = p, qfun = qgbeta, p1 = p1, lower.tail = lower.tail, log.p = log.p,
         shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# random generation
roigbeta <- function(n, shape0, shape1, shape2, p1)
{
  roifun(n = n, rfun = rgbeta, p1 = p1, shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# exposure curve
ecoigbeta <- function(x, shape0, shape1, shape2, p1)
{
  ecoifun(x = x, ecfun = ecgbeta, mfun = mgbeta, p1 = p1, shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# raw moment of a given order
moigbeta <- function(order, shape0, shape1, shape2, p1)
{
  moifun(order = order, mfun = mgbeta, p1 = p1, shape1 = shape1, shape2 = shape2, shape0 = shape0)
}

# total loss probability P(X = 1), i.e. p1
tloigbeta <- function(shape0, shape1, shape2, p1)
{
  tloifun(p1 = p1, shape1 = shape1, shape2 = shape2, shape0 = shape0)
}
