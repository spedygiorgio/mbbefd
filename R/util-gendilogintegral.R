# File overview: legacy code, entirely commented out. It computed the integral
# of log(a + b^x) on [0,1] via the dilogarithm (package gsl) for the variance of the
# MBBEFD distribution. Since version 0.8.14 moments are computed without gsl (see
# mmbbefdR() in distr-mbbefdR-1stparam.R). Kept for reference only.
#
# #integral(from=0, to=1, function(x) a+b^x)
# gendilog <- function(a, b, checkparam=TRUE)
# {
#   if(!(a +1 >0 && b > 0 && a*(1-b) >= 0) && checkparam)
#     return(NaN)
#
#   if(a == 0)
#     return(log(b)/2)
#   else if(b == 1)
#     return(log(a+1))
#   else if(a > 0)
#   {
#     res <- (dilog(-1/a) - dilog(-b/a))/log(b)+log(a)
#   }else
#   {
#     res <- (dilog(1+b/a) - dilog(1+1/a))/log(b)
#     res <- res + log(a+b) - log(-a)*log((a+b)/(a+1))/log(b)
#   }
#   res
# }
