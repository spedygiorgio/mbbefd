library(mbbefd)

data(lossalaefull)
x <- lossalaefull$Loss/lossalaefull$Limit


# fitDR(x, "oigbeta", method="mle", control=list(trace=1, REPORT=1))

dlist <- c("oistpareto", "oibeta", "oigbeta")
flist <- lapply(dlist, function(d) {print(d);
  fitDR(x, d, method="mle")})
names(flist) <- dlist

gof1 <- gofstat(flist)

mm <- rbind(KS = gof1$ks, CvM = gof1$cvm, AD = gof1$ad, AIC = gof1$aic, BIC = gof1$bic)
rownames(mm) <- c("Kolmogorov-Smirnov statistic", "Cramer-von Mises statistic", 
                  "Anderson-Darling statistic", "Aikake's Information Criterion", "Bayesian Information Criterion")


cdfcomp(flist, do.points=FALSE, leg=dlist, addlegend = FALSE, fitlwd=1.5, main="Emp./theo. CDFs - insurance")
legend("topleft", col=c("red", "green", "blue"), leg=dlist, lty=1:3, lwd=2, bty="n")


par(mar=c(4,4,2,1))
cdfcomp(flist, do.points=FALSE, leg=dlist, xlogscale = TRUE, addlegend = FALSE, fitlwd=1.5, main="Emp./theo. CDFs - insurance")
legend("topleft", col=c("red", "green", "blue"), leg=dlist, lty=1:3, lwd=2, bty="n")

eccomp(flist, leg=dlist, do.points = FALSE, addlegend = FALSE, lwd=1.5)
legend("topleft", col=c("red", "green", "blue"), leg=dlist, lty=1:3, lwd=2, bty="n")

ppcomp(flist, leg=dlist, fitpch=".", addlegend = FALSE)
legend("bottomright", fill=c("red", "green", "blue", "cyan", "purple"), leg=dlist)


qqcomp(flist, leg=dlist, use.ppoints=TRUE)
