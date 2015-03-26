require(mbbefd)

#check of expectation for oi distribution
n <- 1e4

probs <- c(1/(2:9))

sapply(probs, function(p) 
{
  x <- roiunif(n, p)
c(mean(x), mbbefd:::tmean1(doiunif, p1=p),
  mbbefd:::tmean2(poiunif, p1=p), mbbefd:::tmean3(x, ecoiunif, p1=p))
}
)



sapply(probs, function(p) 
{
  x <- roistpareto(n, a=2, p)
  c(mean(x), mbbefd:::tmean1(doistpareto, a=2, p1=p),
    mbbefd:::tmean2(poistpareto, a=2, p1=p), mbbefd:::tmean3(x, ecoistpareto, a=2, p1=p))
}
)




sapply(probs, function(p) 
{
  x <- roibeta(n, shape1=2, shape2=3, p)
  c(mean(x), mbbefd:::tmean1(doibeta, shape1=2, shape2=3, p1=p),
    mbbefd:::tmean2(poibeta, shape1=2, shape2=3, p1=p), 
    mbbefd:::tmean3(x, ecoibeta, shape1=2, shape2=3, p1=p))
}
)

