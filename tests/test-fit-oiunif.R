library(mbbefd)
library(fitdistrplus)


#oiunif
n <- 1e3
nboot <- 10
x <- roiunif(n, 1/6)
f1 <- fitdist(x, "oiunif", method="mle", start=list(p1=etl(x)))
b1 <- bootdist(f1, niter=nboot)

summary(b1)

plot(b1)
abline(v=1/6, col="red")

hist(b1$estim[,1])
