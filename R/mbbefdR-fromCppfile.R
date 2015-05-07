

#exposure curve

.G<-function(x, a, b, g)
{
  if(missing(a)){
    # using b, g
    if(identical(g, 1) | identical(b, 0))
      return(x)
    if(identical(b, 1) & g > 1)
      return(
        log(1+(g-1)*x)/log(g)
      )
    if(identical(b*g, 1) & g>1)
      return((1-b^x)/(1-b))
    else
      return(
        log(((g-1)*b + (1-g*b)*b^x)/(1-b))/log(g*b)  
      )
  }else{
    # using a nd b 
    if(identical(a, 0) | identical(b, 1))
      return(x)
    if(a*(1-b)>0)
      return( 
        (log(a+b^x)-log(a+1))/(log(a+b)-log(a+1))
      )
    else
      return((1-b^x) / (1-b))
  }
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

mbbefdExposure<-function(x, a, b, g)
{
  if(x>1||x<0) 
    stop("Error! x should be between 0 and 1")  #check parameters coherence
  if(b<0) 
    stop("b should be greater or equal 0")
  if(missing(a)){
    if(g <= 0) stop("g has to be greater than 0")
    .G(x, g=g,b=b)
  }
  else
    .G(x=x, a=a, b=b)
}

####################################
#classical functions
#distribution function
pmbbefd2<-function(q,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(q>1||q<0) stop("Error! q should lie between 0 and 1")
  out<-1-.Sx(x=q,a=a,b=b)
  return(out)
}

#random generation function (now using the Rcpp version)
#.f4Random<-function(x,a,b) ifelse( ( x>= 1-(a+1)*b/(a+b) ),1,log( (a*(1-x)) / (a+x) ) /log(b))  




#density functio

dmbbefd2<-function(x,a,b,g)
{
  if(missing(a)) a<-g2a(g=g,b=b)
  if(x>1||x<0) stop("Error! x should lie between 0 and 1")
  out <-  -(((a+1)*a*log(b)*b^x)/((a+b^x)^2))*(x<1)+(x==1)*(a+1)*b/(a+b)
  #out <- .dmbbefdC(x=x,a=a,b=b) #TODO: check why it does not work
  return(out)
}


