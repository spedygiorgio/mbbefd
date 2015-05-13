library(mbbefd)

#test of GB1 distribution

#integral of the density
integrate(dgbeta, 0, 1, shape0=1, shape1=3, shape2=3/2)
integrate(dgbeta, 0, 1, shape0=1/2, shape1=3, shape2=3/2)
integrate(dgbeta, 0, 1, shape0=2, shape1=3, shape2=3/2)

z <- 0:10/10
cbind(dbeta(z^(1/pi), 3, 3/2)*1/pi*z^(1/pi-1), dgbeta(z, pi, 3, 3/2))
cbind(pbeta(z^(1/pi), 3, 3/2), pgbeta(z, pi, 3, 3/2))
cbind(qbeta(z, 3, 3/2)^pi, qgbeta(z, pi, 3, 3/2))

#RNG
n <- 1e4
x <- rgbeta(n, shape0=2, shape1=3, shape2=3/2)
y <- rgbeta(n, shape0=pi, shape1=3, shape2=3/2)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), pgbeta(z, shape0=2, shape1=3, shape2=3/2))

cbind(ecdf(y)(z), pgbeta(z, shape0=pi, shape1=3, shape2=3/2))


#mean
c(mean(x), mgbeta(1, shape0=2, shape1=3, shape2=3/2), ecgbeta(1, shape0=2, shape1=3, shape2=3/2))
c(mean(y), mgbeta(1, shape0=pi, shape1=3, shape2=3/2), ecgbeta(1, shape0=pi, shape1=3, shape2=3/2))

#raw moment
for(i in 2:4)
{
  cat("E(X^", i, ")\n", sep="")
  print(c(mean(x^i), mgbeta(i, shape0=2, shape1=3, shape2=3/2)))
  print(c(mean(y^i), mgbeta(i, shape0=pi, shape1=3, shape2=3/2)))
}

#test limited expected value
d <- 1/2
s0 <- 2
s1 <- 3
s2 <- 3/2

mean(pmin(x, d))

f <- function(x, d, shape0) dgbeta(x, shape0=shape0, shape1=3, shape2=3/2)*pmin(x,d)
integrate(f, 0, 1, d=d, shape0=s0)




#test EC
f <- function(x, d, shape0) dgbeta(x, shape0=shape0, shape1=3, shape2=3/2)*pmin(x, d)/mgbeta(1, shape0=shape0, shape1=3, shape2=3/2)
F <- function(d, shape0) sapply(d, function(d) integrate(f, 0, 1, d=d, shape0=shape0)$value)

cbind(eecf(x)(z), ecgbeta(z, shape0=2, shape1=3, shape2=3/2), F(z, shape0=2))

cbind(eecf(y)(z), ecgbeta(z, shape0=pi, shape1=3, shape2=3/2), F(z, shape0=pi))




