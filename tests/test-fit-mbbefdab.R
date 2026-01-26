library(mbbefd)
library(fitdistrplus)
set.seed(123456)

n <- 1e3
nboot <- 1000
nboot <- 10
lossrate <- rmbbefd(n, 1/2, 1/10)


f1 <- fitDR(lossrate, "mbbefd")
summary(f1)

#should be similar
if(FALSE)
{
  f0 <- fitDR(lossrate, "mbbefd", start= list(a=1/4, b=1/4), control=list(trace=1))
  fitDR(lossrate, "mbbefd", start= list(a=1/4, b=1/4), control=list(trace=2))
}



cdfcomp(f1, do.points=FALSE)
qqcomp(f1)
vcov(f1)


# llsurface(plot.min=c(0, 0), plot.max=c(2, 1/2), plot.arg=c("a", "b"), obs=lossrate, distr="mbbefd", nlevels=25)
# points(f1$estimate["a"], f1$estimate["b"], pch="+", col="red")
# points(1/2, 1/10, pch="x", col="black")

b1 <- bootDR(f1, niter=nboot, silent=TRUE)

if(sum(b1$converg == 0) > 2)
{
  plot(b1, enhance=TRUE, trueval=c(1/2, 1/10))
  plot(density(b1))
}

f2 <- fitDR(lossrate, "mbbefd", method="tlmme")
summary(f2)



lossrate <- rmbbefd(n, -1/2, 5)


f1 <- fitDR(lossrate, "mbbefd")
summary(f1)

b1 <- bootDR(f1, niter=nboot, silent=TRUE)
if(sum(b1$converg == 0) > 2)
{
  plot(b1, enhance=TRUE, trueval=c(-1/2, 5))
  plot(density(b1))
}

