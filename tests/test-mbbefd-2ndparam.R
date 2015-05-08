library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

#x <- rMBBEFD(n, 2, 1/2)
x <- rMBBEFD(n, 2, 1/5)
y <- rMBBEFD(n, 3, 2)



#test CDF
z <- 0:8/8
cbind(ecdf(x)(z), pMBBEFD(z, 2, 1/2))

cbind(ecdf(y)(z), pMBBEFD(z, 3, 2))

#test EC
cbind(eecf(x)(z), ecMBBEFD(z, 2, 1/2))

cbind(eecf(y)(z), ecMBBEFD(z, 3, 2))

#test mean
mean(x)
mMBBEFD(1, 2, 1/2)

mean(y)
mMBBEFD(1, 3, 2)

#total loss
etl(x)
tlMBBEFD(2, 1/2)


etl(y)
tlMBBEFD(3, 2)


#test quantile
cbind(quantile(y, probs=0:10/10), qMBBEFD(0:10/10, 3, 2))



z <- seq(0, 1, length=101)
plot(z, pMBBEFD(z, 3, 2), type="l", ylim=c(0, 1-tlMBBEFD(3, 2)))

plot(z, qMBBEFD(z, 3, 2), type="l", xlim=c(0, 1-tlMBBEFD(3, 2)))



#test density

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(x), ylim=c(0,1))
lines(z, dMBBEFD(z, 2, 1/2), col="red")


plot(density(y), ylim=c(0,1.1))
lines(z, dMBBEFD(z, 3, 2), col="red")
