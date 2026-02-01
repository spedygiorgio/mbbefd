library(mbbefd)
library(fitdistrplus)
set.seed(123456)


#oistpareto
n <- 1e3
nboot <- 1000
nboot <- 10
x <- roistpareto(n, 2, 1/6)

f1 <- fitDR(x, "oistpareto", method="mle", control=list(trace=1))
summary(f1)

if(FALSE)
  fitDR(x, "oistpareto", method="mle", control=list(trace=2))
# summary(fitdist(x, "oistpareto", method="mle", start=list(a=1/mean(x), p1=etl(x))))#check

b1 <- bootDR(f1, niter=nboot)
summary(b1)
vcov(b1)
vcov(f1)


if(sum(b1$converg == 0) > 2)
{
  plot(b1, enhance=TRUE, trueval=c(2, 1/6))

  hist(b1$estim[,1])
  hist(b1$estim[,2])
  plot(density(b1))
}

f2 <- fitDR(x, "oistpareto", method="tlmme")
summary(f2)

gofstat(list(f1, f2))
cdfcomp(list(f1, f2), do.points=FALSE)
ppcomp(list(f1, f2))
