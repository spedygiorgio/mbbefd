library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e4

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

