library(mbbefd)
library(fitdistrplus)
set.seed(123456)

#oiunif
n <- 1e3
nboot <- 1000
nboot <- 10
x <- roiunif(n, 1/6)
f1 <- fitDR(x, "oiunif", method="mle", control=list(trace=0))
summary(f1)

if(FALSE)
  fitDR(x, "oiunif", method="mle", control=list(trace=2))

# summary(fitdist(x, "oiunif", method="mle", start=list(p1=1/2))) #check

b1 <- bootDR(f1, niter=nboot)
summary(b1)

if(sum(b1$converg == 0) > 2)
{
  plot(b1)
  abline(v=1/6, col="red")
  
  plot(density(b1))
  
  hist(b1$estim[,1])
  abline(v=1/6, col="red")
}

f2 <- fitDR(x, "oiunif", method="tlmme")


