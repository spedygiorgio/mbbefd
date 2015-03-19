library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e2

x <- rstpareto(n, 2)
y <- rstpareto(n, 1/2)

#test CDF
z <- 0:4/4
ecdf(x)(z)
pstpareto(z, 2)

ecdf(y)(z)
pstpareto(z, 1/2)

#test EC
eecf(x)(z)
ecstpareto(z, 2)

eecf(y)(z)
ecstpareto(z, 1/2)


plot(eecf(x))
v <- seq(0, 1, length=101)
lines(v, ecstpareto(v, 2), lty=3, col="red")


plot(eecf(y))
v <- seq(0, 1, length=101)
lines(v, ecstpareto(v, 1/2), lty=3, col="red")

