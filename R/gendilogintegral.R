
#integral(from=0, to=1, function(x) a+b^x)
gendilog <- function(a, b, checkparam=TRUE)
{
  if(!(a +1 >0 && b > 0 && a*(1-b) >= 0) && checkparam)
    return(NaN)
  
  if(a == 0)
    return(log(b)/2)
  else if(b == 1)
    return(log(a+1))
  require(gsl)
  
  res <- dilog(-1/abs(a)) - dilog(-b/abs(a))
  res <- log(abs(a)) + res/log(b)
  res
}
