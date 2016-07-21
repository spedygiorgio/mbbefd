library(mbbefd)
library(fitdistrplus)
data(smokedfish)

#normed data
listeria <- apply(smokedfish, 1, mean, na.rm=TRUE) / 100
#purety rate
set.seed(1234)
x <- pmin(pmax(jitter(1-listeria, amount=.1), 0.01), 1)

mledist(x[x!=1], "stpareto", start=list(a=2), optim.method="Nelder", control=list(trace=1, REPORT=1), lower=.01)

fitDR(x, "oistpareto", control=list(trace=1, REPORT=1), start=list(a=0.01), optim.method="Nelder")
fitDR(x, "oibeta", control=list(trace=1, REPORT=1))
fitDR(x, "oigbeta", control=list(trace=1, REPORT=1))

dlist <- c("oistpareto", "oibeta", "oigbeta")
flist <- lapply(dlist, function(d) {print(d);
  fitDR(x, d, method="mle")})
names(flist) <- dlist

gof1 <- gofstat(flist)
