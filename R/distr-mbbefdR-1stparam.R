# File overview: d, p, q, r, ec (exposure curve), m (raw moments) and tl (total
# loss) functions of the MBBEFD distribution in its first parametrization (a,b)
# (Bernegger, 1997). For a destruction rate X in [0,1]:
#   exposure curve  G(x) = log((a + b^x) / (a + 1)) / log((a + b) / (a + 1)),
#   survival        1 - F(x) = (a + 1) b^x / (a + b^x)  for 0 <= x < 1,
#   total loss      P(X = 1) = (a + 1) b / (a + b) = 1/g.
# The link with the second parametrization is g = (a + b) / ((a + 1) b) and
# a = (g - 1) b / (1 - g b) (see g2a()).
#
# All functions follow the same structure:
# 1. check the arguments and recycle x, a, b to a common length m;
# 2. start from NaN (invalid parameters), then fill the result case by case:
#    - 'unit indicator' / 'identity function': degenerate parameters, where all
#      the mass is at 1 (total loss only, g = 1);
#    - 'b only': a = +Inf, i.e. gb = 1 in the second parametrization;
#    - 'main case': the regular domain -1 < a < 0 with b > 1, a > 0 with
#      0 < b < 1, or a < -1 with 0 < b < 1 (i.e. b < 1 and gb > 1).
# The user-level names (dmbbefd, ...) are defined at the end of the file.


### R version of d,p,q,r functions MBBEFD(a,b)
#see r functions in distr-mbbefdCpp.R

# Density f(x) = -a (a+1) log(b) b^x / (a + b^x)^2 for 0 <= x < 1, and the
# point mass P(X = 1) = (a+1) b / (a+b) at x = 1 (needed by the likelihood).
dmbbefdR <- function(x, a, b, log=FALSE)
{
  #sanity check
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(x)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(x))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  x <- rep_len(x, length.out=m)
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  id1 <- x == 1
  id01 <- 0 <= x & x < 1 #continuous part on [0,1): f(0) > 0
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  res[idDirac] <- 1 * id1[idDirac] #x == 1
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  res[ididentity] <- 1 * id1[ididentity] #x == 1
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 
  res[idbonly] <- 0 #0
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & id01
  res[idbonly] <- -log(b[idbonly]) * b[idbonly]^x[idbonly]  #-log(b)b^x, x!=1
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & id1
  res[idbonly] <- b[idbonly] #P(X = 1) = b, x==1
  
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  res[idmain] <- 0
  idmain <- idmain & id01
  #-  \ln(b) \frac{a(a+1) b^x }{(a+b^x)^2}
  res[idmain] <- -a[idmain] * (a[idmain]+1) * b[idmain]^x[idmain] * log(b[idmain]) / (a[idmain] + b[idmain]^x[idmain])^2 
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  idmain <- idmain & id1
  #\frac{(a+1)b}{a+b}
  res[idmain] <- (a[idmain]+1) * b[idmain] / (a[idmain]+b[idmain])
  
  if(log)
    res <- log(res)
  
  res
}  

# Distribution function F(q) = 1 - (a+1) b^q / (a + b^q) for 0 < q < 1,
# F(q) = 1 for q >= 1.
pmbbefdR <- function(q, a, b, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(q))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(q)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(q))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  q <- rep_len(q, length.out=m)
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  id0 <- 0 < q 
  id1 <- q < 1
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  res[idDirac] <- 1*(q[idDirac] >= 1)
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  res[ididentity] <- 1 * (q[ididentity] >= 1)
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b)
  res[idbonly] <- 1*(q[idbonly] >= 1) # 1_(x >= 1)
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & id0 & id1
  res[idbonly] <- 1-b[idbonly]^q[idbonly]
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  res[idmain] <- 1*(q[idmain] >= 1) # 1_(x >= 1)
  idmain <- idmain & id0 & id1
  res[idmain] <- 1 - (a[idmain]+1)*b[idmain]^q[idmain]/(a[idmain]+b[idmain]^q[idmain])
  
  if(!lower.tail)
    res <- 1-res
  if(log.p)
    res <- log(res)
  
  res
}  


# Quantile function: for p < 1 - P(X = 1), Q(p) = log((1-p) a / (a + p)) / log(b);
# for larger p the quantile is 1 (total loss).
qmbbefdR <- function(p, a, b, lower.tail = TRUE, log.p = FALSE)
{
  #sanity check
  stopifnot(is.numeric(p))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(p)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(p))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  p <- rep_len(p, length.out=m)
  
  if(!lower.tail)
    p <- 1-p
  if(log.p) 
    p <- exp(p) 
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  idDirac <- idDirac & 0 <= p & p <= 1
  res[idDirac] <- 1
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  ididentity <- ididentity & 0 <= p & p <= 1
  res[ididentity] <- 1
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & 0 <= p & p < 1-b
  res[idbonly] <- log(1 - p[idbonly]) / log(b[idbonly]) #log(1-p)/log(b)
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & p >= 1-b & p <= 1
  res[idbonly] <- 1 #1
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  idmain <- idmain & 0 == p
  res[idmain] <- 0 #0
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  idmain <- idmain & 0 < p & p < 1 - (a+1)*b/(a+b)
  #\frac{\ln\left(\frac{(1-p)a}{a+p}\right)}{\ln(b)}
  res[idmain] <- log( (1 - p[idmain]) * a[idmain] / (a[idmain] + p[idmain]) ) / log( b[idmain] )
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  idmain <- idmain & p >= 1 - (a+1)*b/(a+b) & p <= 1
  res[idmain] <- 1 #1
  
  res
}  


# Vectorized random generation by inversion: Q(U) with U uniform on (0,1).
# Used when the parameters are vectors; the C++ generator is used otherwise
# (see rmbbefd() in distr-mbbefdCpp.R).
rmbbefdR <- function(n, a, b)
{
  #sanity check
  stopifnot(is.numeric(n))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  if(length(n) > 1)
    n <- length(n)
  
  if(min(length(a), length(b), n) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), n)
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  res[idDirac] <- 1
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  res[ididentity] <- 1
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b)
  if(sum(idbonly) > 0)
    res[idbonly] <- qmbbefdR(runif(sum(idbonly)), a[idbonly], b[idbonly])
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(b)
  if(sum(idmain) > 0)
    res[idmain] <- qmbbefdR(runif(sum(idmain)), a[idmain], b[idmain])
  
  res
}


# Exposure curve G(x) = E[min(X,x)] / E[X] (see the file overview).
ecmbbefdR <- function(x, a, b)
{
  #sanity check
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(x)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(x))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  x <- rep_len(x, length.out=m)
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  id0 <- 0 <= x 
  id1 <- x <= 1
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  idDirac <- idDirac & id0 & id1
  res[idDirac] <- 1*(x[idDirac] == 1)
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  ididentity <- ididentity & x > 0 & id1
  res[ididentity] <- x[ididentity]
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & id0 & id1
  res[idbonly] <- (1-b[idbonly]^x[idbonly])/(1-b[idbonly])
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  idmain <- idmain & id0 & id1
  res[idmain] <- log((a[idmain]+b[idmain]^x[idmain])/(a[idmain]+1)) / log((a[idmain]+b[idmain])/(a[idmain]+1))  
  
  res
}

#moment
# Raw moments E[X^k]. The mean has a closed form; other orders use
# E[X^k] = integral on [0,1] of P(X^k > t) dt = integral of S(t^(1/k)) dt,
# with S the survival function, computed by numerical integration.
mmbbefdR <- function(order, a, b)
{
  #sanity check
  stopifnot(is.numeric(order))
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b), length(order)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b), length(order))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  order <- rep_len(order, length.out=m)
  
  res <- rep(NaN, m)
  
  # degenerate case: X = 1, so E[X^k] = 1
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  res[idDirac] <- 1^order[idDirac]
  
  # degenerate case: X = 1
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  res[ididentity] <- 1^order[ididentity]
  
  # a = +Inf: mean (b - 1)/log(b); other orders by integration of S(x) = b^x
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & order == 1
  res[idbonly] <- (b[idbonly] - 1) / log(b[idbonly])
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b) & order != 1
  if(sum(idbonly) > 0)
  {
    surv1 <- function(x, b, k)
      b^(x^(1/k))
    mom1 <- function(b, k)
    {
      res <- try(integrate(surv1, b, k, lower=0, upper = 1))
      if(inherits(res, "try-error"))
        return(NaN)
      res$value
    }
    res[idbonly] <- sapply(1:sum(idbonly), function(i) 
      mom1(k=order[idbonly][i], b=b[idbonly][i]))
  }
  # regular domain: mean (a+1) log((a+b)/(a+1)) / log(b); other orders by integration
  #main case
  idmain <- (-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)
  idmain <- idmain & is.finite(a) & is.finite(b) & order == 1
  res[idmain] <- log( (a[idmain]+b[idmain]) / (a[idmain]+1) ) / log(b[idmain]) * (a[idmain]+1)
  idmain <- (-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)
  idmain <- idmain & is.finite(a) & is.finite(b) & order != 1
  if(sum(idmain) > 0)
  {
    surv2 <- function(x, a, b, k)
      (a+1)*b^(x^(1/k))/(a+b^(x^(1/k)))
    mom2 <- function(a, b, k)
    {
      res <- try(integrate(surv2, a, b, k, lower = 0, upper = 1))
      if(inherits(res, "try-error"))
        return(NaN)
      res$value
    }
    res[idmain] <- sapply(1:sum(idmain), function(i) 
      mom2(k=order[idmain][i], a=a[idmain][i], b=b[idmain][i]))
  }
  
  res
}

#total loss
# Total loss probability P(X = 1) = (a+1) b / (a+b).
tlmbbefdR <- function(a, b)
{
  #sanity check
  stopifnot(is.numeric(a))
  stopifnot(is.numeric(b))
  
  if(min(length(a), length(b)) <= 0)
    return(numeric(0))
  m <- max(length(a), length(b))
  a <- rep_len(a, length.out=m)
  b <- rep_len(b, length.out=m)
  
  #default to NaN: b=+infty; a=+infy && b > 0 && b < 1;
  #main domain: -1 < a < 0 & b > 1, or a > 0 & 0 < b < 1, or a < -1 & 0 < b < 1 (i.e. b < 1 and gb > 1)
  res <- rep(NaN, m)
  
  # degenerate case: total loss only
  #unit indicator - 
  idDirac <- (b == 0 & a > 0) | (a == -1 & b > 1)
  res[idDirac] <- 1
  
  # degenerate case: total loss only (b = 1 or a = 0 gives g = 1)
  #identity function
  ididentity <- (b == 1 & a > -1) | (b == 1 & a < -1) | (b == Inf & a >-1 & a < 0) | (a == 0 & b > 0)
  res[ididentity] <- 1
  
  # a = +Inf (gb = 1): F(x) = 1 - b^x, P(X = 1) = b
  #b only
  idbonly <- a == Inf & b > 0 & b < 1 & is.finite(b)
  res[idbonly] <- b[idbonly]
  
  # regular domain
  #main case
  idmain <- ((-1 < a & a < 0 & b > 1) | (0 < a & 0 < b & b < 1) | (a < -1 & 0 < b & b < 1)) & is.finite(a) & is.finite(b)
  res[idmain] <- (a[idmain] + 1) * b[idmain] / (a[idmain] + b[idmain])
  
  res
}



### d,p,q,ec,m,tl functions MBBEFD(a,b)

# User-level functions (exported): aliases of the R implementations above.
dmbbefd <- dmbbefdR

pmbbefd <- pmbbefdR

qmbbefd <- qmbbefdR

ecmbbefd <- ecmbbefdR

mmbbefd <- mmbbefdR

tlmbbefd <- tlmbbefdR

