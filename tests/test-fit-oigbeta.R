library(mbbefd)
library(fitdistrplus)


#oiunif
n <- 1e2
nboot <- 10
set.seed(12345)
x <- roigbeta(n, 1, 1, 1, 1/6)
x <- rgbeta(n, 1, 1, 1)
m <- mean(x)
v <- (n - 1)/n*var(x)
aux <- m*(1-m)/v - 1
start <- list(shape0=1, shape1=1, shape2=1)


fnobj <- function(par, fix.arg, obs, ddistnam) 
{
  arglist <- c(list(obs), as.list(par), as.list(fix.arg))
  #print(arglist)
  -sum(log(do.call(ddistnam, arglist ) ) )
}

fobj2 <- function(x,y, fix.arg, obs, ddistnam) 
{
  namlist <- names(formals(ddistnam))
  #print(namlist)
  namlist <- namlist[!namlist %in% c("x", "log", names(fix.arg))]
  #print(namlist)
  if(length(namlist) != 2)
    stop("namlist should be a two-element vector.")
  par <- c(x,y)
  names(par) <- namlist
  #print(par)
  fnobj(par, fix.arg, obs, ddistnam)
}


lnL1 <- function(theta, fix.arg)
{
  sapply(theta, function(theta) -fnobj(theta, obs=x, fix.arg=fix.arg, ddistnam=dgbeta))
}
lnL2 <- function(theta, fix.arg)
{
  res <- matrix(nrow=NROW(theta), ncol=NROW(theta))
  for(i in 1:NROW(theta))
      res[i, ] <- sapply(theta[,2], function(y) -fobj2(theta[i,1], y, obs=x, fix.arg=fix.arg, ddistnam=dgbeta))
  res
}

z <- seq(0, 10, by=.05)
fz <- lnL1(z, list(shape1=1, shape2=1))
plot(z, fz, type="l")

fz <- lnL1(z, list(shape1=1, shape0=1))
plot(z, fz, type="l")

fz <- lnL1(z, list(shape0=1, shape2=1))
plot(z, fz, type="l")


z <- seq(0.1, 5, length=51)
fz <- lnL2(cbind(shape0=z,shape1=z), list(shape2=1))
persp(z, z, fz, ticktype="detailed", theta=-120, phi=40)
contour(z, z, fz, levels=c(pretty(range(fz)), -10, -5))

fz <- lnL2(cbind(shape0=z,shape1=z), list(shape1=1))
persp(z, z, fz, ticktype="detailed", theta=-120, phi=40)
contour(z, z, fz, levels=c(pretty(range(fz)), -10, -5))

fz <- lnL2(cbind(shape0=z,shape1=z), list(shape0=1))
persp(z, z, fz, ticktype="detailed", theta=-120, phi=40)
contour(z, z, fz, levels=c(pretty(range(fz)), -10, -5))

#f1 <- fitdist(x, "oigbeta", method="mle", control=list(trace=1, REPORT=1), lower=c(0, 0, 0, 0), 
#              upper=c(Inf, Inf, Inf, 1), optim.method="L-BFGS-B",
#              start=c(start, p1=etl(x)))
f2 <- fitdist(x[x < 1], "gbeta", method="mle", control=list(trace=1, REPORT=1), lower=0.1, 
              upper=Inf, optim.method="L-BFGS-B", start=start)

cdfcomp(f2)

#b1 <- bootdist(f1, niter=nboot)
b2 <- bootdist(f2, niter=nboot)

#summary(b1)
summary(b2)

#plot(b1)
plot(b2)

mbbefd:::pairs2(b2$estim, c(1/2, 3, 2))

hist(b2$estim[,1])
hist(b2$estim[,2])
hist(b2$estim[,3])
