library(mbbefd)

#test of MBBEFD(g,b) distribution
n <- 1e4
set.seed(567)

#### D0 ####

g <- 2 ; b <- 3
obs <- rMBBEFD(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))



#### D1 ####

g <- 2 ; b <- 1/3
obs <- rMBBEFD(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))


#### D2 ####

g <- 2 ; b <- 0
obs <- rMBBEFD(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))


#### D3 ####

g <- 2 ; b <- 1
obs <- mbbefd:::rMBBEFDR(n, g, b)
obs <- rMBBEFD(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))


#### D4 ####

g <- 2 ; b <- Inf
obs <- mbbefd:::rMBBEFDR(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))

#### D5 ####

g <- Inf ; b <- 1/2
obs <- mbbefd:::rMBBEFDR(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))

#### D6 ####

g <- 1 ; b <- 1/2
obs <- mbbefd:::rMBBEFDR(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))

#### D7 ####

g <- 2 ; b <- 1/2
obs <- mbbefd:::rMBBEFDR(n, g, b)

c(tlMBBEFD(g, b), mean(obs==1))

plot(ecdf(obs))
curve(pMBBEFD(x, g, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dMBBEFD(z, g, b), col="blue")

curve(quantile(obs, prob=x))
curve(qMBBEFD(x, g, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecMBBEFD(x, g, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mMBBEFD(k, g, b), sapply(k, function(k) mean(obs^k)))