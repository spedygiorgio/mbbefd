library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4
set.seed(567)

#### D0 ####

a <- -1/2 ; b <- 3
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))


#### D1 ####

a <- 2 ; b <- 1/5
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))



#### D2 ####

a <- 2 ; b <- 0
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))



#### D3 ####

a <- 2 ; b <- 1
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))


#### D4 ####

a <- -1/2 ; b <- Inf
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))


#### D5 ####

a <- -1 ; b <- 3
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))


#### D6 ####

a <- 0 ; b <- 3
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))



#### D7 ####

a <- Inf ; b <- 1/3
obs <- rmbbefd(n, a, b)

c(tlmbbefd(a, b), mean(obs==1))

plot(ecdf(obs))
curve(pmbbefd(x, a, b), add=TRUE, col="blue")

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(obs))
lines(z, dmbbefd(z, a, b), col="blue")

curve(quantile(obs, prob=x))
curve(qmbbefd(x, a, b), add=TRUE, col="blue")

plot(eecf(obs), do.points = FALSE)
curve(ecmbbefd(x, a, b), add=TRUE, col="blue")

k <- 1:10
cbind(k, mmbbefd(k, a, b), sapply(k, function(k) mean(obs^k)))
