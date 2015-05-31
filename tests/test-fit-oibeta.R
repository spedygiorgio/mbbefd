library(mbbefd)
library(fitdistrplus)


#oiunif
n <- 1e3
nboot <- 1000
nboot <- 10
set.seed(12345)
x <- roibeta(n, 3, 2, 1/6)
m <- mean(x)
v <- (n - 1)/n*var(x)
aux <- m*(1-m)/v - 1
start <- list(shape1=m*aux, shape2=(1-m)*aux)

f1 <- fitdist(x, "oibeta", method="mle", control=list(trace=0, REPORT=1), 
			lower=c(0, 0, 0), upper=c(Inf, Inf, 1), optim.method="L-BFGS-B",
            start=c(start, p1=etl(x)))
f2 <- fitdist(x[x < 1], "beta", method="mle", control=list(trace=0, REPORT=1), lower=c(0, 0), 
              upper=c(Inf, Inf), optim.method="L-BFGS-B",
              start=start)

b1 <- bootdist(f1, niter=nboot)
b2 <- bootdist(f2, niter=nboot)

summary(b1)
summary(b2)

plot(b1$estim)

mbbefd:::pairs2(b1$estim, c(3, 2, 1/6), main="one-inflated beta")
mbbefd:::pairs2(b2$estim, c(3, 2))



plot(b2)

hist(b1$estim[,1])
hist(b1$estim[,2])
hist(b1$estim[,3])
