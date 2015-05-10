#include <Rcpp.h>
#include <math.h>  
using namespace Rcpp;


//' Get a parameter known g and b
//' 
//' \code{g2a} returns the a parameter known g and b
//' 
//' @param g the g parameter
//' @param b the b parameter
//' 
//' @return a real value
//' 
//' @examples
//' 
//' g2a(10,2)
//' 
//' @export
// [[Rcpp::export]]
double g2a(double g, double b) {
  double out;
  out=((g-1)*b)/(1-g*b);
  return out;
}

// [[Rcpp::export(.dmbbefdC)]]
double dmbbefdC(double x, double a, double b) {
  double out;
  if (x<1)
    out= -( (a+1)*a*log(b)*pow(b,x) )/( pow( (a+pow(b, x) ),2) );
  else 
    out=((a+1)*b/(a+b));
  return out;
}

//' inverse CDF function
//' 
//' \code{f4Sampler} returns the x known the probability level x and distribution parameters a and b
//' 
//' @param x: the probability
//' @param a: parameter of the mbbefd density function
//' @param b: parameter of the mbbefd density function
//' 
//' @return a real value
//' 
//' @example
//' 
//' f4Sampler(x=.2, a=.2, b=.05)


// [[Rcpp::export(.f4Sampler)]]
double f4Sampler(double x, double a, double b) {
  double out;
  if (x >= 1-(a+1)*b/(a+b))
    out=1;
  else
    out=log( (a*(1-x)) / (a+x) ) /log(b);
  return out;
}


//' random number generation - 1st param
//' 
//' \code{rmbbefdC2} generates random variates distribution parameters a and b
//' 
//' @param n: the number of random variates
//' @param a: first shape parameter
//' @param b: second shape parameter
//' 
//' @return a vector of real values
//' 
//' @example
//' 
//' rmbbefdC2(n=10, a=.2, b=.05)

// [[Rcpp::export(.rmbbefdC2)]]
NumericVector rmbbefdC2(int n, double a, double b) {
  
  NumericVector out(n);
  double u, pab;
  
  if(a == 0 || b == 1) //Dirac
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 0)
        out[i] = 1.0;
      else
        out[i] = 0.0;
    }
  }else if(!isfinite(a))
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 1-b)
        out[i] = 1.0;
      else
        out[i] = log(1-u)/log(b);
    }
  }else
  {
    pab = (a+1)*b/(a+b);
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 1-pab)
        out[i] = 1.0;
      else
        out[i] = log((1-u)*a/(a+u))/log(b);
    }
  }
  
  return out;
}



//' random number generation - 2nd param
//' 
//' \code{rMBBEFDC2} generates random variates distribution parameters g and b
//' 
//' @param n: the number of random variates
//' @param g: first shape parameter
//' @param b: second shape parameter
//' 
//' @return a vector of real values
//' 
//' @example
//' 
//' rMBBEFDC2(n=10, g=2, b=.05)

// [[Rcpp::export(.rMBBEFDC2)]]
NumericVector rMBBEFDC2(int n, double g, double b) {
  
  NumericVector out(n);
  double u;
  
  if(g == 1 || b == 0) //Dirac
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 0)
        out[i] = 1.0;
      else
        out[i] = 0.0;
    }
  }else if(g == 1/b && b < 1) //bg=1
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 1-b)
        out[i] = 1.0;
      else
        out[i] = log(1-u)/log(b);
    }
  }else if(g > 1 && b == 1) //b=1
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 1-1/g)
        out[i] = 1.0;
      else
        out[i] = u/((1-u)*(g-1));
    }
  }else
  {
    for(int i=0; i<n; ++i)
    {
      u = unif_rand();
      if(u > 1-1/g)
        out[i] = 1.0;
      else
        out[i] = 1-log((g*b-1)/(g-1) + (1-b)/((1-u)*(g-1)))/log(b);
    }
  }
  
  return out;
}