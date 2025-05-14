library(mbbefd)

#test of MBBEFD(a,b) distribution
n <- 1e5

# length(b) > 1
a <- 2
b <- c(1/2:10, Inf)
z <- 1/3
cbind(ecmbbefd(z, a, b), log((a+b^z)/(a+1))/log((a+b)/(a+1))) 
cbind(pmbbefd(z, a, b), 1 - (a+1)*b^z/(a+b^z)) 
cbind(dmbbefd(z, a, b), -a * (a+1) * b^z * log(b) / (a + b^z)^2 ) 
cbind(qmbbefd(z, a, b), log((1-z)*a/(a+z))/log(b) ) 
cbind(tlmbbefd(a, b), (a+1)*b/(a+b)) 

#mbbefd:::rmbbefdR(3, a, b)
#rmbbefd(3, a, b)

# length(a) > 1
a <- c(2:10, Inf)
b <- 1/2
z <- 1/3
cbind(ecmbbefd(z, a, b), log((a+b^z)/(a+1))/log((a+b)/(a+1))) 
cbind(pmbbefd(z, a, b), 1 - (a+1)*b^z/(a+b^z)) 
cbind(dmbbefd(z, a, b), -a * (a+1) * b^z * log(b) / (a + b^z)^2 ) 
cbind(qmbbefd(z, a, b), log((1-z)*a/(a+z))/log(b) ) 
cbind(tlmbbefd(a, b), (a+1)*b/(a+b)) 

# length(x) > 1
a <- 2
b <- 1/2
z <- seq(-1/2, 3/2, length=21)
cbind(z, ecmbbefd(z, a, b), log((a+b^z)/(a+1))/log((a+b)/(a+1))) 
cbind(z, pmbbefd(z, a, b), 1 - (a+1)*b^z/(a+b^z)) 
cbind(z, dmbbefd(z, a, b), -a * (a+1) * b^z * log(b) / (a + b^z)^2 ) 
cbind(z, qmbbefd(z, a, b), log((1-z)*a/(a+z))/log(b) ) 

# length(a,b,x) > 1
a <- 2:4
b <- 1/2:5
z <- 1/2:13
a2 <- rep_len(a, length(z))
b2 <- rep_len(b, length(z))
cbind(ecmbbefd(z, a, b), log((a2+b2^z)/(a2+1))/log((a2+b2)/(a2+1))) 
cbind(pmbbefd(z, a, b), 1 - (a2+1)*b2^z/(a2+b2^z)) 
cbind(dmbbefd(z, a, b), -a2 * (a2+1) * b2^z * log(b2) / (a2 + b2^z)^2 ) 
cbind(qmbbefd(z, a, b), log((1-z)*a/(a+z))/log(b) ) 

# length(a,b,x) = 0
a <- 2:4
b <- numeric(0)
z <- 1/2:10
ecmbbefd(z, a, b)
pmbbefd(z, a, b)
dmbbefd(z, a, b)
qmbbefd(z, a, b)

# a, b, x character
a <- "a"
b <- 1/2
z <- 1/2:10
try(ecmbbefd(z, a, b))
try(pmbbefd(z, a, b))
try(dmbbefd(z, a, b))


