# File overview: user-level d, p, q, r, ec (exposure curve), m (raw moment) and
# tl (total loss) functions of the one-inflated shifted truncated Pareto
# distribution: with probability p1 the loss is a total loss (X = 1), otherwise X
# follows the shifted truncated Pareto distribution on [0,1) (see
# distr-shftrncPareto.R). Every function is a thin wrapper around the generic
# one-inflated functions defined in distr-1infl-distrib.R.
# one-inflated shifted truncated Pareto distribution

# density (point mass p1 at x = 1)
doistpareto <- function(x, a, p1, log = FALSE)
{
  doifun(x = x, dfun = dstpareto, p1 = p1, log = log, a = a)
}

# distribution function
poistpareto <- function(q, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  poifun(q = q, pfun = pstpareto, p1 = p1, lower.tail = lower.tail, log.p = log.p, a = a)
}

# quantile function
qoistpareto <- function(p, a, p1, lower.tail = TRUE, log.p = FALSE)
{
  qoifun(p = p, qfun = qstpareto, p1 = p1, lower.tail = lower.tail, log.p = log.p, a = a)
}

# random generation
roistpareto <- function(n, a, p1)
{
  roifun(n = n, rfun = rstpareto, p1 = p1, a = a)
}

# exposure curve
ecoistpareto <- function(x, a, p1)
{
  ecoifun(x = x, ecfun = ecstpareto, mfun = mstpareto, p1 = p1, a = a)
}

# raw moment of a given order
moistpareto <- function(order, a, p1)
{
  moifun(order = order, mfun = mstpareto, p1 = p1, a = a)
}

# total loss probability P(X = 1), i.e. p1
tloistpareto <- function(a, p1)
{
  tloifun(p1 = p1, a = a)
}
