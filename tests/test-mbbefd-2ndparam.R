library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

x <- rMBBEFDR(n, 2, 1/2)
y <- rMBBEFDR(n, 3, 2)



#test CDF
z <- 0:8/8
cbind(ecdf(x)(z), pMBBEFDR(z, 2, 1/2))

cbind(ecdf(y)(z), pMBBEFDR(z, 3, 2))

#test EC
cbind(eecf(x)(z), ecMBBEFDR(z, 2, 1/2))

cbind(eecf(y)(z), ecMBBEFDR(z, 3, 2))

#test mean
mean(x)
mMBBEFDR(1, 2, 1/2)

mean(y)
mMBBEFDR(1, 3, 2)

#total loss
etl(x)
tlMBBEFDR(2, 1/2)


etl(y)
tlMBBEFDR(3, 2)


#test quantile
cbind(quantile(y, probs=0:10/10), qMBBEFDR(0:10/10, 3, 2))



z <- seq(0, 1, length=101)
plot(z, pMBBEFDR(z, 3, 2), type="l", ylim=c(0, 1-tlMBBEFDR(3, 2)))

plot(z, qMBBEFDR(z, 3, 2), type="l", xlim=c(0, 1-tlMBBEFDR(3, 2)))



#test density

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(x), ylim=c(0,1))
lines(z, dMBBEFDR(z, 2, 1/2), col="red")


plot(density(y), ylim=c(0,1.1))
lines(z, dMBBEFDR(z, 3, 2), col="red")
