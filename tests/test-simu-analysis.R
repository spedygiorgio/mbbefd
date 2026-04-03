
library(mbbefd)
set.seed(12345)

nsimu <- 500
nsimu <- 5


#### oigbeta ####

trueval <- list(shape0=3, shape1=4, shape2=3, p1=0.25)

myMLE <- function(n, trueval)
{
  x <- do.call("roigbeta", c(n=n, trueval))
  f1 <- try(fitDR(x, "oigbeta"))
  if(inherits(f1, "fitDR"))
    coef(f1)
  else
    rep(NA, length(trueval))
}


essai200 <- t(replicate(nsimu, myMLE(200, trueval)))
summary(essai200)

#### mbbefd ####


myMLE2 <- function(n, trueval)
{
  x <- do.call("rmbbefd", c(n=n, trueval))
  f1 <- try(fitDR(x, "mbbefd"))
  if(inherits(f1, "fitDR"))
    coef(f1)
  else
    rep(NA, length(trueval))
}


truevalmbbefd <- list(a=3, b=1/2)



essai200 <- t(replicate(nsimu, myMLE2(200, truevalmbbefd)))
summary(essai200)
