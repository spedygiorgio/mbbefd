library(mbbefd)

#test of 1-infl beta distribution
n <- 1e4

x <- roibeta(n, 2, 2, p1=1/2)
y <- roibeta(n, 4/3, 2/3, p1=1/3)

#test pDF
doibeta(1/2, 2, 2, 1/1:10)

#test CDF
z <- 0:10/10
cbind(ecdf(x)(z), poibeta(z, 2, 2, 1/2))

cbind(ecdf(y)(z), poibeta(z, 4/3, 2/3, 1/3))

poibeta(1/2, 2, 2, 1/1:10)

#test quantile
curve(quantile(y, prob=x), n=501, xlab="p")
curve(qoibeta(x, 4/3, 2/3, 1/3), n = 501, add=TRUE, col="green")

#test exposure curve
plot(eecf(y), do.points = FALSE)
curve(ecoibeta(x, 4/3, 2/3, 1/3), n = 501, add=TRUE, col="green")


#total loss
c(etl(x), tloibeta(2, 2, 1/2))
c(etl(y), tloibeta(4/3, 2/3, 1/3))

tloibeta(2, 2, 1/2:10)

#mean
c(mean(x), moibeta(1, 2, 2, 1/2))
c(mean(y), moibeta(1, 4/3, 2/3, 1/3))


#test EC
cbind(eecf(x)(z), ecoibeta(z, 2, 2, 1/2))
cbind(eecf(y)(z), ecoibeta(z, 4/3, 2/3, 1/3))



#plots
n <- 1e2
x <- roibeta(n, 2, 2, p1=1/2)
y <- roibeta(n, 4/3, 2/3, p1=1/3)


plot(eecf(x), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoibeta(v, 2, 2, 1/2), lty=3, col="red")


plot(eecf(y), do.points=FALSE)
v <- seq(0, 1, length=101)
lines(v, ecoibeta(v, 4/3, 2/3, 1/3), lty=3, col="red")

