library(mbbefd)
library(fitdistrplus)


#oigbeta
n <- 1e3
nboot <- 100
nboot <- 10
set.seed(12345)
x <- rgbeta(n, 2, 2, 5/2)

#init value
T.emp <- mean(x/mean(x)*log(x/mean(x)))
psi <- digamma
T.theo  <- function(shape0, shape1, shape2)
{
  EX <- beta(shape1+shape2, 1/shape0) / beta(shape1, 1/shape0)
  1/shape0*(psi(shape1+1/shape0)-psi(shape1+shape2+ 1/shape0)) - log(EX)
}

s0 <- optimize(function(x) (T.emp - T.theo(x, 2, 5/2))^2, lower=0.01, upper=20)$minimum
initpar <- c(list(shape0=s0), as.list(fitdist(x^s0, "beta")$estimate))
par <- list(shape0=2, shape1=2, shape2=5/2)

#test all methods

ctr <- list(trace=0, REPORT=1, maxit=1000)
reslist <- NULL
for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta", optim.method=meth, 
                control=ctr, start=initpar))[3]
  reslist <- c(reslist, list(nograd))
}
for(type in 2:3) #CG with PR or BS updates
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta", optim.method="CG", 
                control=c(ctr, type=type), start=initpar))[3] 
  reslist <- c(reslist, list(nograd)) 
}
fullname <- c("BFGS", "NM", paste("CG", c("FR", "PR", "BS")))
names(reslist) <- fullname

dgbeta2 <- function(x, shape0, shape1, shape2, log=FALSE)
  dgbeta(x, exp(shape0), exp(shape1), exp(shape2), log=log)

initpar2 <- lapply(initpar, log)




for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta2", optim.method="BFGS", 
                          control=ctr, start=initpar2))[3]
  nograd$estimate <- exp(nograd$estimate)
  reslist <- c(reslist, list(nograd))
}
for(type in 2:3) #CG with PR or BS updates
{
  nograd$time <- system.time(nograd <- mledist(x, dist="gbeta2", optim.method="CG", 
                        control=c(ctr, type=type), start=initpar2))[3] 
  nograd$estimate <- exp(nograd$estimate)
  reslist <- c(reslist, list(nograd)) 
}
names(reslist)[(length(fullname)+1):length(reslist)] <- paste("exp.", fullname)

getval <- function(x)
  c(x$estimate, loglik=x$loglik, x$counts, x$time)

resNM <- sapply(reslist[grep("NM", names(reslist))], getval)
resCG <- sapply(reslist[grep("CG", names(reslist))], getval)
resBFGS <- sapply(reslist[grep("BFGS", names(reslist))], getval)  
rownames(resNM) <- rownames(resCG) <- rownames(resBFGS) <- c(paste("fitted", c("shape0", "shape1", "shape2")), "fitted loglik", "func. eval. nb.", "grad. eval. nb.", "time (sec)")


fitdistrplus:::loglikelihood(as.list(nograd$estimate), fix.arg=NULL, obs=x, ddistnam="dgbeta")
nograd$loglik
fitdistrplus:::loglikelihood(as.list(log(nograd$estimate)), fix.arg=NULL, obs=x, ddistnam="dgbeta2")



#empirical check
lnL <- function(par) 
  fitdistrplus:::loglikelihood(as.list(par), fix.arg=NULL, obs=x, ddistnam="dgbeta")

ctr <- list(trace=1, REPORT=1, maxit=2)
dgbeta3 <- function(x, shape0, shape1, shape2, log=FALSE)
{
  cat("x", shape0, shape1, shape2, "\n")
  cat("loglik", lnL(c(shape0, shape1, shape2)), "\n")
  dgbeta(x, shape0, shape1, shape2, log=log)
}

mledist(x, dist="gbeta3", optim.method="BFGS", control=ctr, start=initpar)

fnobj <- function(par, fix.arg, obs, ddistnam){
  -sum(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg), log=TRUE) ) )
}
fnobj(unlist(initpar), NULL, x, "dgbeta")

head(cbind(do.call("dgbeta", c(list(x), initpar, NULL, log=TRUE)), 
log(do.call("dgbeta", c(list(x), initpar, NULL) ) )
))
colSums(cbind(do.call("dgbeta", c(list(x), initpar, NULL, log=TRUE)), 
              log(do.call("dgbeta", c(list(x), initpar, NULL) ) )
))

head(cbind(do.call("dbeta", c(list(x), initpar[2:3], NULL, log=TRUE)), 
           log(do.call("dbeta", c(list(x), initpar[2:3], NULL) ) )
))


lnL(unlist(initpar))





#stochastic algo

library(rgenoud)
mygenoud <- function(fn, par, ...) 
{
  res <- genoud(fn, starting.values=par, ...)        
  c(res, convergence=0)       
}
resgenoud$time <- system.time(resgenoud <- mledist(x, "gbeta", start=initpar, custom.optim= mygenoud, 
            nvars=3, Domains=cbind(c(0,0,0), c(100, 100, 100)), boundary.enforcement=1, 
            hessian=TRUE, print.level=0))[3]

getval(resgenoud)

f1 <- reslist[["BFGS"]]

#test 


f1 <- fitdist(x, "gbeta", method="mle", start=par) # , control=list(trace=3, REPORT=1)) 
summary(f1)
cdfcomp(f1, do.points=FALSE, ylogscale = TRUE)
lines(0:100/100, pgbeta(0:100/100, 2, 2, 5/2), col="green")

denscomp(f1)
lines(0:100/100, dgbeta(0:100/100, 2, 2, 5/2), col="green")


par(mfrow=c(1,3))
llsurface(plot.min=c(0.1, 0.1), plot.max=c(5, 4), nlevels=20,
                         plot.arg=c("shape1", "shape2"), fix.arg=as.list(f1$estimate[1]),
                         plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape1"], f1$estimate["shape2"], pch="+", col="red")
points(2, 5/2, pch="x", col="green")
llsurface(plot.min=c(0.1, 0.1), plot.max=c(6, 4), nlevels=20,
       plot.arg=c("shape0", "shape2"), fix.arg=as.list(f1$estimate[2]),
       plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape0"], f1$estimate["shape2"], pch="+", col="red")
points(2, 5/2, pch="x", col="green")
llsurface(plot.min=c(0.1, 0.1), plot.max=c(5, 6), nlevels=20,
       plot.arg=c("shape1", "shape0"), fix.arg=as.list(f1$estimate[3]),
       plot.np=50, obs=x, distr="gbeta", plot.type="contour")
points(f1$estimate["shape1"], f1$estimate["shape0"], pch="+", col="red")
points(2, 2, pch="x", col="green")


par(mfrow=c(1,3))
llcurve(plot.min=0.1, plot.max=5, plot.arg="shape0", fix.arg=as.list(f1$estimate[-1]), plot.np=50, 
        obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(2, f1$estimate["shape0"]), col=c("green", "red"))
llcurve(plot.min=0.1, plot.max=4, plot.arg="shape1", fix.arg=as.list(f1$estimate[-2]), plot.np=50, 
                      obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(2, f1$estimate["shape1"]), col=c("green", "red"))
llcurve(plot.min=0.1, plot.max=4, plot.arg="shape2", fix.arg=as.list(f1$estimate[-3]), plot.np=50, 
                      obs=x, distr="gbeta", enhance=FALSE)
abline(v=c(5/2, f1$estimate["shape2"]), col=c("green", "red"))


b1 <- bootdist(f1, niter=nboot, silent=TRUE)
summary(b1)

plot(b1, enhance=TRUE, trueval=c(2, 2, 5/2))

constrOptim2 <- function(par, fn, gr=NULL, ui, ci, ...)
  constrOptim(theta=par, f=fn, grad=gr, ui=ui, ci=ci, ...)
f2 <- fitdist(x, "gbeta", start=list(shape0=3, shape1=3, shape2=3), custom.optim=constrOptim2, 
        ui = diag(3), ci = rep(0, 3))
summary(f2)
