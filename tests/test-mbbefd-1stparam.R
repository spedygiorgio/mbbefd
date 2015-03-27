library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

x <- rmbbefdR(n, 2, 1/2)
y <- rmbbefdR(n, -1/2, 2)



#test CDF
z <- 0:8/8
cbind(ecdf(x)(z), pmbbefdR(z, 2, 1/2))

cbind(ecdf(y)(z), pmbbefdR(z, -1/2, 2))

#test EC
cbind(eecf(x)(z), ecmbbefdR(z, 2, 1/2))

cbind(eecf(y)(z), ecmbbefdR(z, -1/2, 2))

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


#test quantile
cbind(quantile(y, probs=0:10/10), qmbbefdR(0:10/10, -1/2, 2))
qmbbefdR(1/2, -1/2, 2)


z <- seq(0, 1, length=101)
plot(z, pmbbefdR(z, -1/2, 2), type="l", ylim=c(0, 1-tlmbbefdR(-1/2, 2)))

plot(z, qmbbefdR(z, -1/2, 2), type="l", xlim=c(0, 1-tlmbbefdR(-1/2, 2)))



#test density

z <- sort(c(1, seq(-0.1,1.1, length=101)))
plot(density(x), ylim=c(0,1))
lines(z, dmbbefdR(z, 2, 1/2), col="red")


plot(density(y), ylim=c(0,1))
lines(z, dmbbefdR(z, -1/2, 2), col="red")
