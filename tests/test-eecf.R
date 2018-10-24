library(mbbefd)

#test of shifted truncated pareto distribution
n <- 1e1

x <- rstpareto(n, 2)+0.01
y <- rstpareto(n, 2)

#test CDF
z <- 0:4/4
ecdf(x)(z)

#test EC
f <- function(d)
  mean(pmin(x, d))/mean(x)
rval <- Vectorize(f, "d")

cbind(eecf(x)(x), rval(x))
cbind(eecf(x)(z), rval(z))

class(eecf(x))
class(ecdf(x))

print(eecf(x))
print(ecdf(x))

cbind(eecf(x)(sort(x)), 
environment(eecf(x))$"Gx")

print(summary(eecf(x)))
print(summary(ecdf(x)))

plot(eecf(x))
plot(ecdf(x))


plot(eecf(x))
plot(eecf(y), add=TRUE, col="red")
lines(eecf(y[1:5]), col="green")

?plot.eecf
?summary.eecf
