library(mbbefd)
library(fitdistrplus)


#oiunif
n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(12345)
x <- roistpareto(n, 2, 1/6)
f1 <- fitdist(x, "oistpareto", method="mle", start=list(a=1/mean(x), p1=etl(x)))
b1 <- bootdist(f1, niter=nboot)

summary(b1)

plot(b1)
abline(v=2, h=1/6, col="red")

hist(b1$estim[,1])
hist(b1$estim[,2])


mbbefd:::pairs2(b1$estim, c(2, 1/6), main="one-inflated ST Pareto")
