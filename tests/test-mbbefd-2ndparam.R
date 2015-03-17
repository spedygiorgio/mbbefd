library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e4

x <- rMBBEFDR(n, 2, 1/2)
y <- rMBBEFDR(n, 4/3, 3)

#test CDF
z <- 0:8/8
ecdf(x)(z)
pMBBEFDR(z, 2, 1/2)

ecdf(y)(z)
pMBBEFDR(z, 4/3, 3) #issue!

#test EC
eecf(x)(z)
ecMBBEFDR(z, 2, 1/2)

eecf(y)(z)
ecMBBEFDR(z, 4/3, 3)

#test mean
mean(x)
mMBBEFDR(1, 2, 1/2)

mean(y)
mMBBEFDR(1, 4/3, 3)

#total loss
etl(x)
tlMBBEFDR(2, 1/2)


etl(y)
tlMBBEFDR(4/3, 3)


