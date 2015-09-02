library(mbbefd)
 
n <- 1e3
nboot <- 100
nboot <- 10
set.seed(123456)
lossrate <- rmbbefd(n, 1/2, 1/10)


f1 <- fitDR(lossrate, "mbbefd")

summary(f1)
cdfcomp(f1, do.points=FALSE)
qqcomp(f1)


b1 <- bootDR(f1, niter=nboot, silent=TRUE)

plot(b1, enhance=TRUE, trueval=c(1/2, 1/10))




f2 <- fitDR(lossrate, "mbbefd", method="tlmme")

