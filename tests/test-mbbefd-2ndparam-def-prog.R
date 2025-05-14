library(mbbefd)

#test of MBBEFD(g,b) distribution
n <- 1e5

# length(b) > 1
g <- 3/2
b <- c(1:6/3, Inf)
z <- 1/4
cbind(g, b, "gb"=g*b, ecMBBEFD(z, g, b), log( ((g-1)*b+(1-g*b)*b^z)/(1-b) ) / log(g*b) )
cbind(g, b, "gb"=g*b, pMBBEFD(z, g, b), 1-(1-b)*b^z/((g-1)*b+(1-g*b)*b^z) )
cbind(g, b, "gb"=g*b, dMBBEFD(z, g, b), -(1-b)*log(b)*(g-1)*b^(1+z)/((g-1)*b+(1-g*b)*b^z)^2 )
cbind(g, b, "gb"=g*b, tlMBBEFD(g, b), 1/g)

# length(x) > 1
g <- 3/2
b <- 1/3
z <- seq(-1/2, 3/2, length=21)
cbind(z, ecMBBEFD(z, g, b), log( ((g-1)*b+(1-g*b)*b^z)/(1-b) ) / log(g*b) )
cbind(z, pMBBEFD(z, g, b), 1-(1-b)*b^z/((g-1)*b+(1-g*b)*b^z) )
cbind(z, dMBBEFD(z, g, b), -(1-b)*log(b)*(g-1)*b^(1+z)/((g-1)*b+(1-g*b)*b^z)^2 )
cbind(z, qMBBEFD(z, g, b))

# length(a,b,x) > 1
g <- 2:4
b <- 1/2:5
z <- 1/2:13
g <- rep_len(g, length(z))
b <- rep_len(b, length(z))

cbind(z, ecMBBEFD(z, g, b), log( ((g-1)*b+(1-g*b)*b^z)/(1-b) ) / log(g*b) )
cbind(z, pMBBEFD(z, g, b), 1-(1-b)*b^z/((g-1)*b+(1-g*b)*b^z) )
cbind(z, dMBBEFD(z, g, b), -(1-b)*log(b)*(g-1)*b^(1+z)/((g-1)*b+(1-g*b)*b^z)^2 )
cbind(g, b, "gb"=g*b, tlMBBEFD(g, b), 1/g)
