library(mbbefd)

data(asiacomrisk)
x <- pmin(asiacomrisk$FGU/asiacomrisk$TIV, 1)
x <- x[!is.na(x)]

plot(ecdf(x))
plot(eecf(x))
etl(x)

fitDR(x, "oibeta", control=list(trace=TRUE))
#fitDR(x, "oigbeta", control=list(trace=TRUE))


dlist <- c("oistpareto", "oibeta", "oigbeta")
dlist <- c("oiunif", "oistpareto", "oibeta")
flist <- lapply(dlist, function(d) {
  print(d);
  fitDR(x, d, method="mle")})
names(flist) <- dlist


cdfcomp(flist, do.points=FALSE, leg=dlist)

ppcomp(flist, leg=dlist, fitpch=".", addlegend = FALSE)
legend("bottomright", fill=c("red", "green", "blue", "cyan"), leg=dlist)

eccomp(flist, leg=dlist, do.points = FALSE)

qqcomp(flist, leg=dlist, use.ppoints=TRUE)
