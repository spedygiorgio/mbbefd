# File overview: log-likelihood of the destruction rate models and, for MBBEFD,
# its analytical gradient and Hessian. They are used by fitDR() for maximum
# likelihood estimation (gradient-based optimizer) and for the standard errors.
#
# Reminder on the likelihood of a destruction rate X in [0,1]: the distribution
# has a density f on [0,1) and a point mass P(X = 1) at total loss. The density
# functions d<dist>() return f(x) for x < 1 and P(X = 1) for x = 1, so that
# LL = sum(log(d<dist>(x))) is the correct log-likelihood for both parts
# (see e.g. Gatti and Wuthrich, 2024, Section 4).

#log-likelihood function
# log-likelihood: dist is the distribution name (the density is d<dist>),
# theta the vector of parameters in the order of the density arguments
LLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD", "unif", "stpareto", "beta", "gbeta"))
  ddist <- paste0("d", dist)
  sum(log(do.call(ddist, c(list(obs), as.list(theta)) ) ) )
}


#gradient of the log-likelihood function : only valid for parameter domain D1, D2
# Gradient of the log-likelihood with respect to the two parameters.
# The formulas are valid on all the domains used by fitDR() (b > 1; b < 1 with
# gb < 1; b < 1 with gb > 1), as they only differentiate log|.| terms.
# Each g1 / g2 helper returns the derivative of log f(x) for one observation;
# the results are summed over the observations.
grLLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("mbbefd", "MBBEFD")) 
  # MBBEFD(a,b): for x < 1, f(x) = -a (a+1) log(b) b^x / (a + b^x)^2 and
  # for x = 1, P(X = 1) = (a+1) b / (a + b).
  if(dist == "mbbefd")
  {
    # d log f / d a
    g1 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, (b-1)/(a+1)/(a+b), (2*a+1)/(a*(a+1)) - 2/(a+b^x))
    }
    # d log f / d b
    g2 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, a/(b*(a+b)), x/b+1/(b*log(b))-2*b^x*x/(b*(a+b^x)))
    }
    c(sum(sapply(obs, g1, theta=theta)), sum(sapply(obs, g2, theta=theta)))
  # MBBEFD(g,b): for x < 1, f(x) = (g-1)(b-1) log(b) b^(1-x) / ((g-1) b^(1-x) + 1 - gb)^2
  # and for x = 1, P(X = 1) = 1/g.
  }else
  {
    # d log f / d g
    g1 <- function(x, theta)
    {
      g <- theta[1]; b <- theta[2]
      
      if(x != 1)
      {
        denom1 <- (g-1)*b^(1-x)+1-g*b
        res <- 1/(g-1) -2*(b^(1-x) - b)/denom1
      }else
      {
        res <- -1/g
      }
      res
    }
    # d log f / d b
    g2 <- function(x, theta)
    {
      g <- theta[1]; b <- theta[2]
      
      if(x != 1)
      {
        denom1 <- (g-1)*b^(1-x)+1-g*b
        num1 <- 1/(b-1)+1/(b*log(b))+(1-x)/b
        #derivative of -2*log(denom1) w.r.t. b
        num2 <- -2*((g-1)*(1-x)*b^(-x) - g)/denom1
        res <- num1+num2
      }else
      {
        res <- 0
      }
      res
    }
    c(sum(sapply(obs, g1, theta=theta)), sum(sapply(obs, g2, theta=theta)))
  }
}

#Hessian of the log-likelihood function
# Hessian of the log-likelihood (matrix of second derivatives), only for
# MBBEFD(a,b). h11 = d2/da2, h21 = d2/dadb, h22 = d2/db2 for one observation,
# summed over the observations. It is used by fitDR() to compute the standard
# errors of MBBEFD(a,b) fits; it is checked against finite differences in the tests.
heLLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("mbbefd", "MBBEFD")) 
  if(dist == "mbbefd")
  {
    h11 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2-1/(a+1)^2, 2/(a+b^x)^2-1/a^2-1/(a+1)^2)
    }
    h21 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2, 2*x*b^(x-1)/(a+b^x)^2)
    }
    h22 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, 1/(a+b)^2-1/b^2, 
             x/b^2-(log(b)+1)/(b^2*log(b)^2)-2*a*x/(b^2*(a+b^x))-2*a*x^2*b^x/(b^2*(a+b^x)^2))
    }
    rbind(c(sum(sapply(obs, h11, theta=theta)), sum(sapply(obs, h21, theta=theta))),
          c(sum(sapply(obs, h21, theta=theta)), sum(sapply(obs, h22, theta=theta))))
  }else
  {
    stop("not yet implemented.")
  }
}
