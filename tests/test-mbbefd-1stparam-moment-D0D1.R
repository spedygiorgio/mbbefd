

library(mbbefd)

#proposition 5 - moment d'ordre k
f <- function(x, k,  a, b)
  pmbbefd(x^(1/k), a, b, lower.tail = FALSE)

integrate(f, k=1, a=2, b=1/2, lower=0, upper = 1)
mmbbefd(1, 2, 1/2)

integrate(f, k=2, a=2, b=1/2, lower=0, upper = 1)
mmbbefd(2, 2, 1/2)



