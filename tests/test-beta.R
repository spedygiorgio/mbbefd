library(mbbefd)

#test of beta distribution
n <- 1e4

x <- roibeta(n, 2, 2, p1=0)
y <- roibeta(n, 4/3, 2/3, p1=0)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), poibeta(z, 2, 2, 0))

cbind(ecdf(y)(z), poibeta(z, 4/3, 2/3, 0))

#total loss
c(etl(x), tloibeta(2, 2, 0))
c(etl(y), tloibeta(4/3, 2/3, 0))

#mean
c(mean(x), moibeta(1, 2, 2, 0))
c(mean(y), moibeta(1, 4/3, 2/3, 0))


#test EC
cbind(eecf(x)(z), ecoibeta(z, 2, 2, 0))
cbind(eecf(y)(z), ecoibeta(z, 4/3, 2/3, 0))



#plots
n <- 1e5
n <- 1e2
x <- roibeta(n, 2, 2, p1=0)
y <- roibeta(n, 4/3, 2/3, p1=0)


plot(eecf(x), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoibeta(v, 2, 2, 0), lty=3, col="red")


plot(eecf(y), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoibeta(v, 4/3, 2/3, 0), lty=3, col="red")

