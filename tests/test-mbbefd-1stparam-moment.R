

library(mbbefd)

f <- function(x, k,  a, b)
  pmbbefd(x^(1/k), a, b, lower.tail = FALSE)

integrate(f, k=1, a=2, b=1/2, lower=0, upper = 1)
mmbbefd(1, 2, 1/2)

integrate(function(x, b) b^(x), b=1/2, lower=0, upper = 1)
mmbbefd(1, Inf, 1/2)

integrate(f, k=2, a=2, b=1/2, lower=0, upper = 1)
mmbbefd(2, 2, 1/2)

integrate(function(x, b) b^(sqrt(x)), b=1/2, lower=0, upper = 1)
mmbbefd(2, Inf, 1/2)

#integrate(f, k=pi, a=2, b=1/2, lower=0, upper = 1)
#mmbbefd(pi, 2, 1/2)

integrate(function(x, b, k) b^(x^(1/pi)), k=pi, b=1/2, lower=0, upper = 1)
#mmbbefd(pi, Inf, 1/2)

