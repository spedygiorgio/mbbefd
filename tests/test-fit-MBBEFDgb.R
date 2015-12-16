library(mbbefd)
 
n <- 1e3
nboot <- 100
nboot <- 10
set.seed(123456)
x <- rMBBEFD(n, 8, 1/4)

f1 <- fitDR(x, "MBBEFD", control=list(trace=1, REPORT=1))

summary(f1)
cdfcomp(f1, do.points=FALSE)
qqcomp(f1)



llsurface(plot.min=c(1, 0), plot.max=c(10, 1/2), plot.arg=c("g", "b"), obs=x, distr="MBBEFD", nlevels=25)
points(f1$estimate["a"], f1$estimate["b"], pch="+", col="red")
points(8, 1/4, pch="x", col="black")


#b1 <- bootDR(f1, niter=nboot, silent=TRUE)
#plot(b1, enhance=TRUE, trueval=c(2, 1/10))

