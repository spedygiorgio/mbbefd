library(mbbefd)
library(fitdistrplus)
data(smokedfish)

#normed data
listeria <- apply(smokedfish, 1, mean, na.rm=TRUE) / 100
#purety rate
set.seed(1234)
x <- pmin(pmax(jitter(1-listeria, amount=.2), 0.01), 1)

#mledist(x[x!=1], "stpareto", start=list(a=2), optim.method="Nelder", control=list(trace=1, REPORT=1), lower=.01)

#, control=list(trace=1, REPORT=1)
flist <- list(
  fitDR(x, "oistpareto", start=list(a=0.01), optim.method="Nelder"),
  fitDR(x, "oibeta"), fitDR(x, "oigbeta"))
names(flist) <- dlist <- c("oistpareto", "oibeta", "oigbeta")

gof1 <- gofstat(flist)

mm <- rbind(KS = gof1$ks, CvM = gof1$cvm, AD = gof1$ad, AIC = gof1$aic, BIC = gof1$bic)
rownames(mm) <- c("Kolmogorov-Smirnov statistic", "Cramer-von Mises statistic", 
                  "Anderson-Darling statistic", "Aikake's Information Criterion", "Bayesian Information Criterion")


cdfcomp(flist, do.points=FALSE, leg=dlist, addlegend = FALSE, lwd=1.5, main="Emp./theo. CDFs - ecotoxicology")
legend("topleft", col=c("red", "green", "blue"), leg=dlist, lty=1:3, lwd=2, bty="n")

#log scale
par(mar=c(4,4,2,1))
cdfcomp(flist, do.points=FALSE, leg=dlist, xlogscale = TRUE, addlegend = FALSE, lwd=1.5, main="Emp./theo. CDFs - ecotoxicology")
legend("topleft", col=c("red", "green", "blue"), leg=dlist, lty=1:3, lwd=2, bty="n")
