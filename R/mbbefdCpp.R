#' @useDynLib mbbefd
#' @importFrom Rcpp sourceCpp




#random generation function (now using the Rcpp version)
#.f4Random<-function(x,a,b) ifelse( ( x>= 1-(a+1)*b/(a+b) ),1,log( (a*(1-x)) / (a+x) ) /log(b))  


#inverse distribution funtcion (quantile)
qmbbefdC<-function(p,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(p>1||p<0) stop("Error! p should lie between 0 and 1")
  #out<- .f4Random(x = p,a=a,b=b) #this s the R native version
  out<- .f4Sampler(x = p,a=a,b=b) #now using the Rcpp version instead
  #p > 1 - prob total loss would make to exceed 1. So it will return one.
  out<-ifelse(out>1, 1,out)
  return(out)
}

#random generation
rmbbefdC <- function(n, a, b, g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  out<-numeric(n)
  u<-runif(n=n,min=0,max=1)
  out<-sapply(u,.f4Sampler,a=a,b=b) #using the Rcpp function
  return(out)
}

rMBBEFDC <- function(n, g, b)
{
  a<-g2a(g=g,b=b)
  out<-numeric(n)
  u<-runif(n=n,min=0,max=1)
  out<-sapply(u,.f4Sampler,a=a,b=b) #using the Rcpp function
  return(out)
}



### r function MBBEFD(a,b)
rmbbefd <- rmbbefdC

### r function MBBEFD(g,b)
rMBBEFD <- rMBBEFDC

