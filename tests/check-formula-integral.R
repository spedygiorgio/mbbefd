
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