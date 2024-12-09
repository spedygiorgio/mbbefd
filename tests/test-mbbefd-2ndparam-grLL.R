library(mbbefd)

g1 <- function(g, b, x)
{
  denom1 <- (g-1)*b^(1-x)+1-g*b
  num1 <- (b-1)*log(b)*b^(1-x) -2*(b-1)*(g-1)*log(b)*b^(1-x)*(b^(1-x)-b)/denom1
  num1/denom1^2
}
G1 <- function(x, g, b)
{
  sapply(g, function(y) integrate(g1, lower=1, upper=y, b=b, x=x)$value)
}


G1(1/2, 3:10, 2)
sapply(3:10, function(g) dMBBEFD(1/2, g, 2))


g2 <- function(b, g, x)
{
  denom1 <- (g-1)*b^(1-x)+1-g*b
  num1 <- (g-1)*log(b)*b^(1-x) + (b-1)*(g-1)*b^(-x) 
  num2 <- (b-1)*(g-1)*log(b)*(1-x)*b^(-x) 
  num3 <- - 2*(b-1)*(g-1)*log(b)*b^(1-x) *((g-1)*(1-x)*b^(-x)-g)/denom1
  (num1+num2+num3)/denom1^2
}

g21 <- function(b, g, x)
{
  U <- (g-1)*b^(1-x)+1-g*b
  Uprimex <- -(g-1)*log(b)*b^(1-x)
  T1 <- -(1-b)*(1/b/log(b) - (1-x)/b)
  T2 <- 2*(1-b)/U*( (g-1)*(1-x)*b^(-x) -g )
  -Uprimex/U^2*(1 + T1 + T2)
}



G2 <- function(x, g, b)
{
  sapply(b, function(y) integrate(g2, lower=1, upper=y, g=g, x=x)$value)
}


G21 <- function(x, g, b)
{
  sapply(b, function(y) integrate(g21, lower=1, upper=y, g=g, x=x)$value)
}

G22 <- function(x, g, b)
{
  sapply(b, function(y) try(integrate(g22, lower=1.1, upper=y, g=g, x=x)$value))
}

x0 <- 1/3
cbind(G2(x0, 3, 2:10),
  G21(x0, 3, 2:10),
  G22(x0, 3, 2:10),
  sapply(2:10, function(b) dMBBEFD(x0, 3, b))
)

