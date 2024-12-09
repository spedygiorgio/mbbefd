
#log-likelihood function
LLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("oiunif", "oistpareto", "oibeta", "oigbeta", "mbbefd", "MBBEFD", "unif", "stpareto", "beta", "gbeta"))
  ddist <- paste0("d", dist)
  sum(log(do.call(ddist, c(list(obs), as.list(theta)) ) ) )
}

#gradient of the log-likelihood function
grLLfunc <- function(obs, theta, dist)
{
  dist <- match.arg(dist, c("mbbefd", "MBBEFD")) 
  if(dist == "mbbefd")
  {
    g1 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, (b-1)/(a+1)/(a+b), (2*a+1)/(a*(a+1)) - 2/(a+b^x))
    }
    g2 <- function(x, theta)
    {
      a <- theta[1]; b <- theta[2]
      ifelse(x == 1, a/(b*(a+b)), x/b+1/(b*log(b))-2*b^x*x/(b*(a+b^x)))
    }
    c(sum(sapply(obs, g1, theta=theta)), sum(sapply(obs, g2, theta=theta)))
  }else
  {
    g2 <- function(x, theta)
    {
      g <- theta[1]; b <- theta[2]
      
      if(x != 1)
      {
        denom1 <- (g-1)*b^(1-x)+1-g*b
        num1 <- (g-1)*log(b)*b^(1-x) + (b-1)*(g-1)*b^(-x) 
        num2 <- (b-1)*(g-1)*log(b)*(1-x)*b^(-x) 
        num3 <- - 2*(b-1)*(g-1)*log(b)*b^(1-x) *((g-1)*(1-x)*b^(-x)-g)/denom1
        res <- (num1+num2+num3)/denom1^2
      }else
      {
        res <- 1
      }
      res
    }
    g1 <- function(x, theta)
    {
      g <- theta[1]; b <- theta[2]
      
      denom1 <- (g-1)*b^(1-x)+1-g*b
      num1 <- (b-1)*log(b)*b^(1-x) -2*(b-1)*(g-1)*log(b)*b^(1-x)*(b^(1-x)-b)/denom1
      num1/denom1^2
    }
    c(sum(sapply(obs, g1, theta=theta)), sum(sapply(obs, g2, theta=theta)))
  }
}

#Hessian of the log-likelihood function
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
