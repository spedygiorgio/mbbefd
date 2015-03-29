#' @useDynLib mbbefd
#' @importFrom Rcpp sourceCpp

# 
# g2a<-function(g,b) {
#   out<-((g-1)*b)/(1-g*b)
#   return(out)
# }



#exposure curve

.G<-function(x,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  out<-(log(a+b^x)-log(a+1))/(log(a+b)-log(a+1))
  return(out)
}

#its derivative
dG<-function(x,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  out<-( (log(b)*b^x) /(a+b^x) )/(log(a+b)-log(a+1))
  return(out)
}

#the survival function
.Sx<-function(x,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  out<-dG(x=x,a=a,b=b)/dG(x=0,a=a,b=b)
  return(out)
}


#the function to compute the exposure function

mbbefdExposure<-function(x, a, b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(x>1||x<0) stop("Error! x should be between 0 and 1")  #check parameters coherence
  if(b<0) stop("b should be greater or equal 0")
  out<-.G(x=x, a=a, b=b)
  return(out)
}

####################################
#classical functions
#distirbution function
pmbbefd<-function(q,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(q>1||q<0) stop("Error! q should lie between 0 and 1")
  out<-1-.Sx(x=q,a=a,b=b)
  return(out)
}

#random generation function (now using the Rcpp version)
#.f4Random<-function(x,a,b) ifelse( ( x>= 1-(a+1)*b/(a+b) ),1,log( (a*(1-x)) / (a+x) ) /log(b))  


#inverse distribution funtcion (quantile)
qmbbefd<-function(p,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(p>1||p<0) stop("Error! p should lie between 0 and 1")
  #out<- .f4Random(x = p,a=a,b=b) #this s the R native version
  out<- .f4Sampler(x = p,a=a,b=b) #now using the Rcpp version instead
  #p > 1 - prob total loss would make to exceed 1. So it will return one.
  out<-ifelse(out>1, 1,out)
  return(out)
}

#density functio

dmbbefd<-function(x,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(x>1||x<0) stop("Error! x should lie between 0 and 1")
  out <-  -(((a+1)*a*log(b)*b^x)/((a+b^x)^2))*(x<1)+(x==1)*(a+1)*b/(a+b)
  #out <- .dmbbefdC(x=x,a=a,b=b) #TODO: check why it does not work
  return(out)
}


#min(1, (log( (a*(1-x)) / (a+x) )) /log(b) )
# 
# b=12.648
# g=4.22069

rmbbefd<-function(n,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  out<-numeric(n)
  u<-runif(n=n,min=0,max=1)
  out<-sapply(u,.f4Sampler,a=a,b=b) #using the Rcpp function
  return(out)
}



####Swiss Re curves####

swissRe<-function(c)
{
  out<-numeric(2)
  b <- exp(3.1 - 0.15*c*(1+c))
  g <-exp(c*(0.78 + 0.12*c))
  out<-c(b,g)
  names(out)<-c("b","g")
  return(out)
}