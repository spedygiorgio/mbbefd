library(mbbefd)

data(lossalaefull)
lossalaefull <- na.omit(lossalaefull)

x <- lossalaefull$Loss/lossalaefull$Limit

mean(doigbeta(x, 0.7216387, 0.2962170, 1.5459744, log = TRUE, p1=1/2))

#### manual optim ####

f <- function(theta)
  -mean(doigbeta(x, theta[1], theta[2], theta[3], log = TRUE, p1=theta[4]))

optim(c(0.7216387, 0.2962170, 1.5459744, 1/2), f, lower=0.001, method="L-BFGS-B")$par
optim(c(0.7216387, 0.2962170, 1.5459744, 1/2), f)$par

optim(c(0.7216387, 0.2962170, 1.5459744, 1/2), f, lower=c(0.001, 1, 0.001, 0.001), method="L-BFGS-B")$par
optim(c(0.7216387, 0.2962170, 1.5459744, 1/2), f, lower=c(0.001, 1, 0.001, 0.001))$par


#### fitdist ####

fitdist(x[x < 1], "gbeta", method="mle", 
        start=list(shape0=0.7216387, shape1=0.2962170, shape2=1.5459744),
        control=list(trace=0))
fitdist(x[x < 1], "gbeta", method="mle", lower=0, optim.method="L-BFGS-B",
        start=list(shape0=0.7216387, shape1=0.2962170, shape2=1.5459744),
        control=list(trace=0))

fitdist(x[x < 1], "gbeta", method="mle", lower=0, optim.method="L-BFGS-B",
        start=list(shape0=0.02130868, shape1=39.27109204, shape2=3.16571652),
        control=list(trace=0))

fitdist(x, "oigbeta", method="mle", lower=0, optim.method="L-BFGS-B",
        start=list(shape0=0.7216387, shape1=0.2962170, shape2=1.5459744, p1=5/100),
        control=list(trace=0))


fitDR(x, "oigbeta", method="mle", control=list(trace=2, REPORT=10))


