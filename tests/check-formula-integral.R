library(mbbefd)

#Formula 1
f <- function(x) log(x)/x
integrate(f, 1/4, 1)
-log(1/4)^2/2

#Formula 2 -> D2
library(gsl)
f <- function(x) log(1-x)/x
g <- 3
b <- 1/4
integrate(f, b*(g-1)/(g*b-1), (g-1)/(g*b-1))
-dilog((g-1)/(g*b-1)) + dilog(b*(g-1)/(g*b-1))

#Formula 2 -> D1
g <- 3
b <- 4
integrate(f, b*(g-1)/(g*b-1), (g-1)/(g*b-1))
-dilog((g-1)/(g*b-1)) + dilog(b*(g-1)/(g*b-1))

#formula 2 -> D2

f <- function(x) log(1+x*(g-1)/(1-g*b))/x
g <- 3
b <- 1/4
integrate(f, b, 1)
-dilog((g-1)/(g*b-1)) + dilog(b*(g-1)/(g*b-1)) 

#formula 2 -> D2

f <- function(x) log(1-g*b+x*(g-1))/x
g <- 3
b <- 1/4
integrate(f, b, 1)
-dilog((g-1)/(g*b-1)) + dilog(b*(g-1)/(g*b-1)) - log(1-g*b)*log(b)



#Formula 3
f <- function(x, g, b)
  log(x)/x/( (g-1)*x+1-g*b )

I <- function(g, b)
{
  temp <- dilog(b*(g-1)/(g*b-1)) - dilog((g-1)/(g*b-1)) 
  (log(b)*log(abs(1-b)) - log(b)^2/2 +temp - log(abs(1-g*b))*log(b) )/(1-g*b)
}

b <- 1/4
integrate(f, b, 1, g=3, b=b)$value
I(g=3, b=b)

b <- 4
integrate(f, b, 1, g=3, b=b)$value
I(g=3, b=b)


#Formula second order moment D3
f <- function(x)
    1/(1+(g-1)*sqrt(x))

I <- function(g)
    2/(g-1)-2*log(g)/(g-1)^2

g <- 3.5
integrate(f, 0, 1)
I(3.5)


#formula of regularity conditions - Lemma B1

f <- function(x)
  (a+b^x)^m 

m <- 1; a <- 3; b <- 1/2
integrate(f, 0, 1)
m <- 2; a <- 3; b <- 1/2
integrate(f, 0, 1)
m <- 3; a <- 3; b <- 1/2
integrate(f, 0, 1)

I <- function(m, a, b)
{
  if(m == 0)
    return(1)
  if(m > 0)
  {
    k <- 1:m
    return(a^m+ sum(choose(m,k)*a^(m-k)/log(b)*(b^k/k - 1/k)))
  }else
  {
    m <- -m 
    k <- 1:(m-1)
    res <- 1/a^{m} - {log({a+b}/{a+1})}/{a^{m}*log(b)}
    if(m >= 2)
    res <- res + sum({-1}/{a^k*log(b)*(-m+k)}*({1}/{(a+b)^{m-k}} - {1}/{(a+1)^{m-k}}))
    return(res)
  }
}
I(1, a, b)
I(2, a, b)
I(3, a, b)

g <- function(x)
  (a+b^x)^(-m) 
m <- 1; a <- 3; b <- 1/2
integrate(g, 0, 1)
I(-1, a, b)
m <- 2; a <- 3; b <- 1/2
integrate(g, 0, 1)
I(-2, a, b)
m <- 3; a <- 3; b <- 1/2
integrate(g, 0, 1)
I(-3, a, b)



J <- function(m, a, b)
{
  L2ab <- mbbefd:::gendilog(a,b)
  if(m == 0)
    return(1)
  if(m == 1)
    return(1/(2*a)-(log(a+b) - L2ab)/(a*log(b)))
  if(m == 2)
    return(1/a^2/2-log(a+b)/a^2/log(b)+log((a+b)/(a+1))/a^2/log(b)^2 + L2ab/a^2/log(b)-b/log(b)/a^2/(a+b))
  k <- 1:m
  J(m-1, a, b)/a - 1/log(b)/a/(-m+1)/(a+b)^{m-1} + I(-m+1, a, b)/log(b)/a/(-m+1)
}

g <- function(x)
  x/(a+b^x)^(m) 

m <- 1; a <- 3; b <- 1/2
integrate(g, 0, 1)
J(1, a, b)
m <- 2; a <- 3; b <- 1/2
integrate(g, 0, 1)
J(2, a, b)
m <- 3; a <- 3; b <- 1/2
integrate(g, 0, 1)
J(3, a, b)
m <- 4; a <- 3; b <- 1/2
integrate(g, 0, 1)
J(4, a, b)


#formula of regularity conditions - Lemma B2
f <- function(x)
  b^x*log(b)/((a+b^x)^m)

I <- function(m, a, b)
{
  if(m == 1)
    return(log((a+b)/(a+1)))
  1/(-m+1)*(1/(a+b)^{m-1} - 1/(a+1)^{m-1})
}

m <- 1; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(1, a, b)

m <- 2; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(2, a, b)


#formula of regularity conditions - Lemma B3


f <- function(x)
  x*b^x*log(b)/(a+b^x)^m



I <- function(m, a, b)
{
  if(m == 1)
    return(log(a+b)-mbbefd:::gendilog(a,b))
  res <- b/a/(a+b) - log((a+b)/(a+1))/a/log(b)
  if(m>2)
  {
    l <- 1:(m-2)
    res <- res - sum(1/a^l/log(b)/(m-1)/(-m+1+l)*(1/(a+b)^(m-1-l) - 1/(a+1)^(m-1-l)) ) 
  }
  res  
}
  
m <- 1; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)

m <- 2; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)

m <- 3; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)

m <- 4; a <- 3; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)


#formula of regularity conditions - Lemma B4

f <- function(x)
  x^2*b^x*log(b)/(a+b^x)^m

I <- function(m, a, b)
{
  L2ab <- mbbefd:::gendilog(a,b)
  
  if(m == 2)
    res <- -1/(a+b)+1/a+2/a/log(b)*(-log(a+b)+L2ab)
  else if(m == 3)
    res <- -1/2/(a+b)^2+1/a^2/2-log(a+b)/a^2/log(b)+log((a+b)/(a+1))/a^2/log(b)^2 + L2ab/a^2/log(b)-b/log(b)/a^2/(a+b)
  else
    res <- NA
  res  
}

m <- 2; a <- 3.1; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)

m <- 3; a <- 3.1; b <- 1/2
integrate(f, 0, 1)
I(m, a,b)
