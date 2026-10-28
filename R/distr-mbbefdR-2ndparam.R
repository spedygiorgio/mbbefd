# File overview: d, p, q, r, ec (exposure curve), m (raw moments) and tl (total
# loss) functions of the MBBEFD distribution in its second parametrization (g,b)
# (Bernegger, 1997), with g >= 1 and b >= 0. For a destruction rate X in [0,1]:
#   P(X = 1) = 1/g (probability of a total loss),
#   general case (b != 1, gb != 1):
#     exposure curve G(x) = log(((g-1) b + (1 - gb) b^x) / (1 - b)) / log(gb),
#     1 - F(x) = (1 - b) / ((g-1) b^(1-x) + 1 - gb) for 0 <= x < 1;
#   special case b = 1 ('g only'): G(x) = log(1 + (g-1) x) / log(g);
#   special case gb = 1 ('b only'): G(x) = (1 - b^x) / (1 - b).
# The general case includes the three domains used by fitDR(): b > 1;
# b < 1 with gb < 1; b < 1 with gb > 1 (see Gatti and Wuthrich, 2024).
#
# All functions follow the same structure:
# 1. check the arguments and recycle x, g, b to a common length m;
# 2. start from NaN (invalid parameters), then fill the result case by case:
#    'unit indicator' / 'identity function' (degenerate: total loss only, g = 1,
#    or b = 0 / b = +Inf), 'b only' (gb = 1), 'g only' (b = 1), 'main case'.
# The user-level names (dMBBEFD, ...) are defined at the end of the file.


### R version of d,p,q,r functions MBBEFD(g,b)
#see r functions in distr-mbbefdCpp.R

# Density for 0 <= x < 1 and point mass P(X = 1) = 1/g at x = 1
# (needed by the likelihood).
dMBBEFDR <- function(x, g, b, log = FALSE)
{
  #sanity check
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b), length(x)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(x))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)
  x <- rep_len(x, length.out = m)

  #default to NaN
  res <- rep(NaN, m)
  id1 <- x == 1
  id01 <- 0 <= x & x < 1 #continuous part on [0,1): f(0) > 0

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1 & is.finite(b)
  res[idDirac] <- 1 * id1[idDirac] #x == 1

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  res[ididentity] <- 1 * id1[ididentity] #x == 1

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b)
  res[idbonly] <- 0
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & id01
  res[idbonly] <- -log(b[idbonly]) * b[idbonly]^x[idbonly]  #-log(b)b^x, x!=1
  idbonly <- g > 1 & b != 1 & b > 0 & b * g == 1 & is.finite(b) & id1
  res[idbonly] <- b[idbonly] #P(X = 1) = b, x==1

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- g > 1 & b == 1
  res[idgonly] <- 0
  idgonly <- g > 1 & b == 1 & is.finite(g) & id01
  #(g-1)/(1+(g-1)x)^2
  res[idgonly] <- (g[idgonly] - 1) / (1 + (g[idgonly] - 1) * x[idgonly])^2
  idgonly <- g > 1 & b == 1 & is.finite(g) & id1
  #1/g
  res[idgonly] <- 1 / g[idgonly]

  # general case b != 1, gb != 1
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b)
  res[idmain] <- 0
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & id01
  # - \frac{(1-b)\ln(b)(g-1)b^{1+x} }{ ( (g-1)b+(1-gb)b^x)^2}
  res[idmain] <- -(1 - b[idmain]) * (g[idmain] - 1) * log(b[idmain]) * b[idmain]^(1 - x[idmain])
  res[idmain] <- res[idmain] / ((g[idmain] - 1) * b[idmain]^(1 - x[idmain]) + 1 - g[idmain] * b[idmain])^2
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & id1
  #1/g
  res[idmain] <- 1 / g[idmain]

  if (log)
    res <- log(res)
  res
}

# Distribution function; F(q) = 1 for q >= 1.
pMBBEFDR <- function(q, g, b, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(q))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b), length(q)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(q))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)
  q <- rep_len(q, length.out = m)

  #default to NaN
  res <- rep(NaN, m)
  id0 <- 0 < q
  id1 <- q < 1

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1
  res[idDirac] <- 1 * (q[idDirac] >= 1) # 1_(x >= 1)

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  res[ididentity] <- 1 * (q[ididentity] >= 1) # 1_(x >= 1)

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1  & is.finite(b)
  res[idbonly] <- 1 * (q[idbonly] >= 1) # 1_(x >= 1)
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & id0 & id1
  res[idbonly] <- 1 - b[idbonly]^q[idbonly] #1-b^x

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- g > 1 & b == 1 & is.finite(g)
  res[idgonly] <- 1 * (q[idgonly] >= 1) # 1_(x >= 1)
  idgonly <- g > 1 & b == 1 & is.finite(g) & id0 & id1
  #1-1/(1+(g-1)x)
  res[idgonly] <- 1 - 1 / (1 + (g[idgonly] - 1) * q[idgonly])

  # general case b != 1, gb != 1
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b)
  res[idmain] <- 1 * (q[idmain] >= 1) # 1_(x >= 1)
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & id0 & id1
  #1-(1-b)/((g-1)b^(1-x) + 1-gb)
  res[idmain] <- 1 - (1 - b[idmain]) / ((g[idmain] - 1) * b[idmain]^(1 - q[idmain]) + 1 - g[idmain] * b[idmain])

  if (!lower.tail)
    res <- 1 - res
  if (log.p)
    res <- log(res)

  res
}


# Quantile function: inverse of F below 1 - P(X = 1), and 1 (total loss) above.
qMBBEFDR <- function(p, g, b, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b), length(p)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(p))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)
  p <- rep_len(p, length.out = m)

  #default to NaN
  res <- rep(NaN, m)

  if (!lower.tail)
    p <- 1 - p
  if (log.p)
    p <- exp(p)

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1
  idDirac <- idDirac & 0 <= p & p <= 1
  res[idDirac] <- 1

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  ididentity <- ididentity & 0 <= p & p <= 1
  res[ididentity] <- 1

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & 0 <= p & p < 1 - b
  res[idbonly] <- log(1 - p[idbonly]) / log(b[idbonly]) #log(1-p)/log(b)
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & p >= 1 - b & p <= 1
  res[idbonly] <- 1 #1

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- g > 1 & b == 1 & is.finite(g) & 0 <= p & p < 1 - 1 / g
  res[idgonly] <- p[idgonly] / (1 - p[idgonly]) / (g[idgonly] - 1) #p/(1-p)/(g-1)
  idgonly <- g > 1 & b == 1 & is.finite(g) & p >= 1 - 1 / g & p <= 1
  res[idgonly] <- 1 #1

  # general case b != 1, gb != 1
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & 0 == p
  #0
  res[idmain] <- 0
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & 0 < p & p < 1 - 1 / g
  #(gb-1)/(g-1) + (1-b)/(1-p)/(g-1)
  res[idmain] <- (g[idmain] * b[idmain] - 1) / (g[idmain] - 1) + (1 - b[idmain]) / (1 - p[idmain]) / (g[idmain] - 1)
  #1-log(.)/log(b)
  res[idmain] <- 1 - log(res[idmain]) / log(b[idmain])
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & p >= 1 - 1 / g & p <= 1
  res[idmain] <- 1 #1

  res
}


# Vectorized random generation by inversion: Q(U) with U uniform on (0,1).
# Used when the parameters are vectors; the C++ generator is used otherwise
# (see rMBBEFD() in distr-mbbefdCpp.R).
rMBBEFDR <- function(n, g, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))
  if (length(n) > 1)
    n <- length(n)

  if (min(length(g), length(b), n) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), n)
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)

  #default
  res <- rep(NaN, m)

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  # (comparisons with NA give NA: %in% TRUE turns them into FALSE, so that
  # missing parameters give NaN)
  idDirac <- (g == Inf & b > 0 & b != 1) %in% TRUE
  res[idDirac] <- 1

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- ((g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)) %in% TRUE
  res[ididentity] <- 1

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- (g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b)) %in% TRUE
  if (sum(idbonly) > 0)
    res[idbonly] <- qMBBEFDR(runif(sum(idbonly)), g[idbonly], b[idbonly])

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- (g > 1 & b == 1 & is.finite(g)) %in% TRUE
  if (sum(idgonly) > 0)
    res[idgonly] <- qMBBEFDR(runif(sum(idgonly)), g[idgonly], b[idgonly])

  # general case b != 1, gb != 1
  #main case
  idmain <- (g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b)) %in% TRUE
  if (sum(idmain) > 0)
    res[idmain] <- qMBBEFDR(runif(sum(idmain)), g[idmain], b[idmain])

  res
}


# Exposure curve G(x) = E[min(X,x)] / E[X] (see the file overview).
ecMBBEFDR <- function(x, g, b)
{
  #sanity check
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b), length(x)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(x))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)
  x <- rep_len(x, length.out = m)

  #default to NaN
  res <- rep(NaN, m)
  id0 <- 0 <= x
  id1 <- x <= 1

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1 & id0 & id1
  res[idDirac] <- 1 * (x[idDirac] > 0) # 1_(0 < x < 1)

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  ididentity <- ididentity & id0 & id1
  res[ididentity] <- x[ididentity] #x

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & id0 & id1
  res[idbonly] <- (1 - b[idbonly]^x[idbonly]) / (1 - b[idbonly]) #(1-b^x)/(1-b)

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- g > 1 & b == 1 & is.finite(g) & id0 & id1
  #log(1+(g-1)x)/log(g)
  res[idgonly] <- log(1 + (g[idgonly] - 1) * x[idgonly]) / log(g[idgonly])

  # general case b != 1, gb != 1
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & id0 & id1
  #(g-1)b+(1-gb)b^x
  res[idmain] <- (g[idmain] - 1) * b[idmain] + (1 - g[idmain] * b[idmain]) * b[idmain]^x[idmain]
  #log(./(1-b))/log(gb)
  res[idmain] <- log(res[idmain] / (1 - b[idmain])) / log(g[idmain] * b[idmain])

  res
}

#moment
# Raw moments E[X^k]. The mean has a closed form in each case; other orders
# use E[X^k] = integral on [0,1] of S(t^(1/k)) dt (S = survival function),
# computed by numerical integration (surv*/mom* helpers below).
mMBBEFDR <- function(order, g, b)
{
  #sanity check
  stopifnot(is.numeric(order))
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b), length(order)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b), length(order))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)
  order <- rep_len(order, length.out = m)

  res <- rep(NaN, m)

  # degenerate case: X = 1, so E[X^k] = 1
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1
  res[idDirac] <- 1^order[idDirac]

  # degenerate case: X = 1
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  res[ididentity] <- 1^order[ididentity]

  # gb = 1: mean (b - 1)/log(b)
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & order == 1
  res[idbonly] <- (b[idbonly] - 1) / log(b[idbonly])
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b) & order != 1
  if (sum(idbonly) > 0)
  {
    surv3 <- function(x, b, k)
      b^(x^(1 / k))
    mom3 <- function(b, k)
    {
      res <- try(integrate(surv3, b, k, lower = 0, upper = 1))
      if (inherits(res, "try-error"))
        return(NaN)
      res$value
    }
    res[idbonly] <- sapply(1:sum(idbonly), function(i)
      mom3(k = order[idbonly][i], b = b[idbonly][i]))
  }

  # b = 1: mean log(g)/(g - 1)
  #g only
  idgonly <- g > 1 & b == 1 & is.finite(g) & order == 1
  res[idgonly] <- log(g[idgonly]) / (g[idgonly] - 1)
  idgonly <- g > 1 & b == 1 & is.finite(g) & order != 1
  if (sum(idgonly) > 0)
  {
    surv4 <- function(x, g, k)
      1 / (1 + (g - 1) * x^(1 / k))
    mom4 <- function(g, k)
    {
      res <- try(integrate(surv4, g, k, lower = 0, upper = 1))
      if (inherits(res, "try-error"))
        return(NaN)
      res$value
    }
    res[idgonly] <- sapply(1:sum(idgonly), function(i)
      mom4(k = order[idgonly][i], g = g[idgonly][i]))
  }

  # general case: mean (b-1) log(gb) / ((gb-1) log(b)) (written equivalently below)
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & order == 1
  res[idmain] <- (1 - b[idmain]) * log(g[idmain] * b[idmain]) / (1 - g[idmain] * b[idmain]) / log(b[idmain])
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b) & order != 1
  if (sum(idmain) > 0)
  {
    surv5 <- function(x, g, b, k)
      (1 - b) / ((g - 1) * b^(1 - x^(1 / k)) + 1 - g * b)
    mom5 <- function(g, b, k)
    {
      res <- try(integrate(surv5, g, b, k, lower = 0, upper = 1))
      if (inherits(res, "try-error"))
        return(NaN)
      res$value
    }
    res[idmain] <- sapply(1:sum(idmain), function(i)
      mom5(k = order[idmain][i], g = g[idmain][i], b = b[idmain][i]))
  }
  res
}


#total loss
# Total loss probability P(X = 1): 1/g in the regular cases.
tlMBBEFDR <- function(g, b)
{
  #sanity check
  stopifnot(is.numeric(g))
  stopifnot(is.numeric(b))

  if (min(length(g), length(b)) <= 0)
    return(numeric(0))
  m <- max(length(g), length(b))
  g <- rep_len(g, length.out = m)
  b <- rep_len(b, length.out = m)

  #default to NaN
  res <- rep(NaN, m)

  # degenerate case: total loss only (g = 1)
  #unit indicator -
  idDirac <- g == Inf & b > 0 & b != 1
  res[idDirac] <- 1

  # degenerate case: total loss only (b = 0, b = +Inf or g = 1)
  #identity function
  ididentity <- (g > 1 & b == 0) | (g > 1 & b == Inf) | (g == 1 & b > 0 & b != 1)
  res[ididentity] <- 1

  # special case gb = 1: F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- g > 1 & b < 1 & b > 0 & b * g == 1 & is.finite(b)
  res[idbonly] <- b[idbonly]

  # special case b = 1: F(x) = 1 - 1/(1 + (g-1) x), P(X = 1) = 1/g
  #g only
  idgonly <- g > 1 & b == 1 & is.finite(g)
  res[idgonly] <- 1 / g[idgonly]

  # general case b != 1, gb != 1
  #main case
  idmain <- g > 1 & b > 0 & b != 1 & b * g != 1 & is.finite(g) & is.finite(b)
  res[idmain] <- 1 / g[idmain]

  res
}




### d,p,q,ec,m,tl functions MBBEFD(g,b)

# User-level functions (exported): aliases of the R implementations above.
dMBBEFD <- dMBBEFDR

pMBBEFD <- pMBBEFDR

qMBBEFD <- qMBBEFDR

ecMBBEFD <- ecMBBEFDR

mMBBEFD <- mMBBEFDR

tlMBBEFD <- tlMBBEFDR
