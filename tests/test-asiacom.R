library(mbbefd)

data(asiacomrisk)
x <- asiacomrisk$DR
x <- x[!is.na(x)]

plot(ecdf(x))
plot(eecf(x))
etl(x)

#test optim method
if(FALSE)
{
  fitDR(x, "oistpareto", control=list(trace=TRUE))  
  
fitDR(x, "oibeta", control=list(trace=TRUE))
fitDR(x, "oibeta", control=list(trace=TRUE), optim.method="L-BFGS-B")

fitDR(x, "oigbeta", control=list(trace=TRUE))
fitDR(x, "oigbeta", control=list(trace=TRUE), optim.method="L-BFGS-B")
fitDR(x, "oigbeta", control=list(trace=TRUE), optim.method="BFGS")
}

dlist <- c("oiunif", "oistpareto", "oibeta", "oigbeta")
flist <- lapply(dlist, function(d) {
  cat("distribution:", d, "\n");
  fitDR(x, d, method="mle", optim.method=ifelse(d=="oigbeta", "BFGS", "default"))})
names(flist) <- dlist


cdfcomp(flist, do.points=FALSE, leg=dlist)

ppcomp(flist, leg=dlist, fitpch=".", addlegend = FALSE)
legend("bottomright", fill=c("red", "green", "blue", "cyan"), leg=dlist)

eccomp(flist, leg=dlist, do.points = FALSE)

qqcomp(flist, leg=dlist, use.ppoints=TRUE)
