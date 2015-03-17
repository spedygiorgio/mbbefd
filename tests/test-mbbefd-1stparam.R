library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

x <- rmbbefdR(n, 2, 1/2)
y <- rmbbefdR(n, -1/2, 2)

#test CDF
z <- 0:8/8
ecdf(x)(z)
pmbbefdR(z, 2, 1/2)

ecdf(y)(z)
pmbbefdR(z, -1/2, 2)

#test EC
eecf(x)(z)
ecmbbefdR(z, 2, 1/2)

eecf(y)(z)
ecmbbefdR(z, -1/2, 2)

#test mean
mean(x)
mmbbefdR(1, 2, 1/2)

mean(y)
mmbbefdR(1, -1/2, 2)

#total loss
etl(x)
tlmbbefdR(2, 1/2)


etl(y)
tlmbbefdR(-1/2, 2)


