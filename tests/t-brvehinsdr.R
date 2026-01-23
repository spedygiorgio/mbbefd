library(mbbefd)

x <- read.csv("~/Desktop/brvehinsdr.csv")[,1]
f1 <- fitDR(x, "oibeta", method="tlmme", control=list(trace=1))
coef(f1)

f1 <- fitDR(x, "oibeta", method="mle", control=list(trace=1))
coef(f1)


f2 <- fitDR(x, "oigbeta", method="tlmme", control=list(trace=1))
coef(f2)



DIFF2 <- function(par, obs) 
{
  PX1 <- do.call(paste0("tl", dist), as.list(par))
  EX <- do.call(paste0("m", dist), as.list(c(order=1, par)))
  if(npar <= 2)
    return( (EX - mean(obs))^2 + (PX1 - etl(obs))^2 )
  
  if(npar >= 3)
    EX2 <- do.call(paste0("m", dist), as.list(c(order=2, par)))
  if(npar >= 4)
    EX3 <- do.call(paste0("m", dist), as.list(c(order=3, par)))
  
  if(npar == 3)
    return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (PX1 - etl(obs))^2 )
  else if(npar == 4)
    return( (EX - mean(obs))^2 + (EX2 - mean(obs^2))^2 + (EX3 - mean(obs^3))^2 + (PX1 - etl(obs))^2 )
  else
    stop("not implemented")
}

dist <- "oibeta"

x00 <- list(shape1=0.245913, shape2=1.383036, p1= 0.248)
npar <- length(x00)

DIFF2(unlist(x00), x)


optim(unlist(x00), DIFF2, obs=x)
optim(c(1, 1, 1/2), DIFF2, obs=x)
optim(unlist(x00), DIFF2, obs=x, lower=0, upper=c(Inf, Inf, 1), method="L-BFGS-B", control=list(trace=1, REPORT=1))



