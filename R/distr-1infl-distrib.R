# File overview: generic functions for one-inflated distributions on [0,1].
# A one-inflated distribution mixes a total loss with probability p1 (X = 1)
# and a continuous distribution on [0,1) with probability 1 - p1:
#   F(x) = (1 - p1) F0(x) + p1 * 1(x >= 1).
# Each function takes the corresponding function of the continuous part
# (dfun, pfun, qfun, rfun, ecfun, mfun) and passes its parameters through "...".
# All functions are vectorized in x and p1 (recycled to a common length), and
# return NaN for invalid p1 (outside [0,1]) or missing values.
#d, p, q, r function for one-inflated distribution


# "Density": (1 - p1) f0(x) for 0 <= x < 1, p1 at x = 1 (the point mass,
# as needed by the likelihood), 0 outside [0,1].
doifun <- function(x, dfun, p1, log=FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(x))
  p1 <- rep_len(p1, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  
  #x == 1
  idmain <- x == 1 & p1 >= 0 & p1 <= 1 & !is.na(x) & !is.na(p1)
  if(!log)
  {
    res[idmain] <- p1[idmain]  
  }else
    res[idmain] <- log(p1[idmain])
  # x not in [0,1]
  idmain <- (x > 1 | x < 0) & p1 >= 0 & p1 <= 1 & !is.na(x) & !is.na(p1)
  if(!log)
  {
    res[idmain] <- 0
  }else
    res[idmain] <- -Inf
  # x in [0,1)
  idmain <- 0 <= x & x < 1 & p1 >= 0 & p1 <= 1 & !is.na(x) & !is.na(p1)
  if(!log)
  {
    res[idmain] <- dfun(x[idmain], log=FALSE, ...) * (1 - p1[idmain])
  }else
    res[idmain] <- dfun(x[idmain], log=TRUE, ...) + log(1 - p1[idmain])
  res
}

# Distribution function (1 - p1) F0(q) + p1 1(q >= 1).
poifun <- function(q, pfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(q)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(q))
  p1 <- rep_len(p1, length.out=m)
  q <- rep_len(q, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1 & !is.nan(q) & !is.na(p1)
  res[idmain] <- pfun(q[idmain], lower.tail = TRUE, log.p = FALSE, ...)*(1 - p1[idmain]) + p1[idmain]*(q[idmain] >= 1)
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}


# Quantile function: for p < 1 - p1 it is the quantile of the continuous part at
# level p / (1 - p1); for p >= 1 - p1 it is 1 (total loss).
qoifun <- function(p, qfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(p)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(p))
  p1 <- rep_len(p1, length.out=m)
  p <- rep_len(p, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 < 1 & p >= 0 & p < 1 - p1 & !is.na(p) & !is.na(p1)
  p[idmain] <- p[idmain]/(1-p1[idmain]) #transformed quantile
  if(!lower.tail)
    p[idmain] <- 1-p[idmain]
  if(log.p) 
    p[idmain] <- exp(p[idmain]) 
  
  res[idmain] <- qfun(p[idmain], lower.tail = TRUE, log.p = FALSE, ...)
  
  idmain <- p1 >= 0 & p1 < 1 & p <= 1 & p >= 1 - p1 & !is.na(p) & !is.na(p1)
  res[idmain] <- 1
  res
}  

# Random generation: draw a Bernoulli(p1) total loss indicator, then draw
# from the continuous part for the observations that are not total losses.
roifun <- function(n, rfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  if(length(n) > 1)
    n <- length(n)
  if(min(length(p1), n) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), n)
  p1 <- rep_len(p1, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1 & !is.na(p1)
  res[idmain] <- rbinom(sum(idmain), 1, p1[idmain])
  
  notone <- res[idmain] != 1 & idmain
  res[notone] <- rfun(sum(notone), ...)
  
  res
}

#exposure curve and moment functions
# Exposure curve of the mixture. With G0 and E0 the exposure curve and mean of
# the continuous part, E[min(X,x)] = (1-p1) E0 G0(x) + p1 x and E[X] = (1-p1) E0 + p1,
# so G(x) = ((1-p1) G0(x) + p1 x / E0) / (1 - p1 + p1 / E0).
ecoifun <- function(x, ecfun, mfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  
  if(min(length(p1), length(x)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(x))
  p1 <- rep_len(p1, length.out=m)
  x <- rep_len(x, length.out=m)
  
  res <- rep(NaN, m)
  
  idmain <- x >= 0 & x <= 1  & p1 >= 0 & p1 <= 1 & !is.na(x) & !is.na(p1)
  
  G0 <- ecfun(x[idmain], ...) #exposure curve
  E0 <- mfun(order=1, ...) #expectation
  
  res[idmain] <- ((1-p1[idmain])*G0[idmain] + p1[idmain]*x[idmain]/E0)/(1-p1[idmain]+p1[idmain]/E0)
  res
}


# moment function
# Raw moment E[X^k] = p1 * 1^k + (1 - p1) E0[X^k].
moifun <- function(order, mfun, p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  stopifnot(is.numeric(order))
  
  if(min(length(p1), length(order)) <= 0)
    return(numeric(0))
  
  m <- max(length(p1), length(order))
  p1 <- rep_len(p1, length.out=m)
  order <- rep_len(order, length.out=m)
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1 & !is.na(order) & !is.na(p1)
  E0 <- mfun(order=order[idmain], ...) #expectation
  res[idmain] <- p1[idmain] + (1-p1[idmain])*E0
  res
}

#total loss function
# Total loss probability P(X = 1) = p1.
tloifun <- function(p1, ...)
{
  #sanity check
  stopifnot(is.numeric(p1))
  if(length(p1) <= 0)
    return(numeric(0))
  
  m <- length(p1)
  res <- rep(NaN, m)
  
  idmain <- p1 >= 0 & p1 <= 1 & !is.na(p1)
  res[idmain] <- p1[idmain]
  res
}
