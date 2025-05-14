

library(mbbefd)

f <- function(x, k,  g, b)
  pMBBEFD(x^(1/k), g, b, lower.tail = FALSE)

integrate(f, k=1, g=2, b=1/3, lower=0, upper = 1)
mMBBEFD(1, 2, 1/3)

integrate(f, k=2, g=2, b=1/3, lower=0, upper = 1)
mMBBEFD(2, 2, 1/3)

integrate(f, k=pi, g=2, b=1/3, lower=0, upper = 1)
mMBBEFD(pi, 2, 1/3)



integrate(function(x, b) b^(x), b=1/2, lower=0, upper = 1)
mMBBEFD(1, 2, 1/2)

integrate(function(x, b) b^(sqrt(x)), b=1/2, lower=0, upper = 1)
mMBBEFD(2, 2, 1/2)

integrate(function(x, b, k) b^(x^(1/pi)), k=pi, b=1/2, lower=0, upper = 1)
mMBBEFD(pi, 2, 1/2)

