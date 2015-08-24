library(mbbefd)
 
n <- 1e3
nboot <- 1000 
set.seed(123456)
lossrate <- rmbbefd(n, 1/2, 1/10)

library(fitdistrplus)


initparmbbefd <- list(list(a=-1/2, b=2), list(a=2, b=1/2))
distlist <- "mbbefd"



#to meet the standard 'fn' argument and specific name arguments, we wrap optimize,
library(alabama)
args(constrOptim.nl)

constrab <- function(x, fix.arg, obs, ddistnam)
{
	x[1]*(1-x[2]) #a*(1-b) >= 0
}

f1 <- fitDR(lossrate, "mbbefd")

summary(f1)
cdfcomp(f1, do.points=FALSE)
qqcomp(f1)

alabama1 <- fitdist(lossrate, distr=distlist, start=initparmbbefd[[1]], 
        custom.optim= constrOptim.nl, hin=constrab, method="mle",
        control.outer=list(trace= FALSE))

mbbefd:::grLLfunc(obs=lossrate, theta=c(0.5548252, 0.1069368), dist="mbbefd")
mbbefd:::heLLfunc(obs=lossrate, theta=c(0.5548252, 0.1069368), dist="mbbefd")
  
alabama2 <- fitdist(lossrate, distr=distlist, start=initparmbbefd[[2]], 
        custom.optim= constrOptim.nl, hin=constrab, method="mle",
        control.outer=list(trace=FALSE))


#solution at the frontier for (a,b) in (-1, 0)x(1, +infty)
c(alabama2$estimate, alabama1$estimate)
c(alabama2$loglik, alabama1$loglik)

b1 <- bootdist(alabama1, niter=nboot)
b2 <- bootdist(alabama2, niter=nboot)

mbbefd:::pairs2(b1$estim, c(1/2, 1/10), main="MBBEFD(a,b)")
mbbefd:::pairs2(b1$estim, main="MBBEFD(a,b)", trueval=NULL)
mbbefd:::pairs2(b2$estim, c(1/2, 1/10), main="MBBEFD(a,b)")
