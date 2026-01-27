library(mbbefd)
library(fitdistrplus)
set.seed(123456)

n <- 1e3
nboot <- 100
nboot <- 10
x <- rMBBEFD(n, 8, 1/4)

system.time(f1 <- fitDR(x, "MBBEFD"))
summary(f1)

#should be similar
if(FALSE)
{
  f0 <- fitDR(x, "MBBEFD", start= list(g=4, b=1/2), control=list(trace=1))
  fitDR(x, "MBBEFD", start= list(g=4, b=1/2), control=list(trace=2))
  fitDR(x, "MBBEFD", start= list(g=4, b=1/2), method="tlmme", control=list(trace=2))
}

cdfcomp(f1, do.points=FALSE)
qqcomp(f1)



# llsurface(plot.min=c(1, 0), plot.max=c(11, 1/2), plot.arg=c("g", "b"), obs=x, distr="MBBEFD", nlevels=25)
# points(f1$estimate["g"], f1$estimate["b"], pch="+", col="red")
# points(8, 1/4, pch="x", col="black")


b1 <- bootDR(f1, niter=nboot, silent=TRUE)
if(sum(b1$converg == 0) > 2)
{
  plot(b1, enhance=TRUE, trueval=c(8, 1/4))
  plot(density(b1))
}


if(FALSE)
{
  x <- rMBBEFD(n, 2, 1/4)
  
  system.time(f1 <- fitDR(x, "MBBEFD"))
  summary(f1)
  
  b1 <- bootDR(f1, niter=nboot, silent=TRUE)
  plot(b1, enhance=TRUE, trueval=c(2, 1/4))
}
