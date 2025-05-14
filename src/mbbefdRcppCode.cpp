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
  out = ( ( g - 1 ) * b) / ( 1 - g * b );
  return out;
}


//' random number generation - 1st param
//' 
//' \code{rmbbefdC} generates random variates distribution parameters a and b
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

// [[Rcpp::export(.rmbbefdC)]]
NumericVector rmbbefdC(int n, double a, double b) {
  
  NumericVector out(n);
  double u, pab;
  if(a == 0 || b == 1 || !R_FINITE(b) || b == 0 || a == -1) //Dirac
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 0)
        out[i] = 1.0;
      else
        out[i] = 0.0;
    }
  }else if(!R_FINITE(a))
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-b)
        out[i] = 1.0;
      else
        out[i] = log(1-u)/log(b);
    }
  }else if(a +1 >0 && b > 0 && a*(1-b) > 0)
  {
    pab = (a+1)*b/(a+b);
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-pab)
        out[i] = 1.0;
      else
        out[i] = log((1-u)*a/(a+u))/log(b);
    }
  }else 
  {
    for(int i=0; i<n; i++)
    {
      out[i] = R_NaN;
    }
  } 
  
  return out;
}



//' random number generation - 2nd param
//' 
//' \code{rMBBEFDC} generates random variates distribution parameters g and b
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

// [[Rcpp::export(.rMBBEFDC)]]
NumericVector rMBBEFDC(int n, double g, double b) {
  
  NumericVector out(n);
  double u;
  
  if(g == 1 || !R_FINITE(g) || b == 0 || !R_FINITE(b)) //Dirac
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 0)
        out[i] = 1.0;
      else
        out[i] = 0.0;
    }
  }else if(b*g == 1 && b < 1) //bg=1
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-b)
        out[i] = 1.0;
      else
        out[i] = log(1-u)/log(b);
    }
  }else if(g > 1 && b == 1) //b=1
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-1/g)
        out[i] = 1.0;
      else
        out[i] = u/((1-u)*(g-1));
    }
  }else if(b > 0 & g > 1)
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-1/g)
        out[i] = 1.0;
      else
        out[i] = 1-log((g*b-1)/(g-1) + (1-b)/((1-u)*(g-1)))/log(b);
    }
  }else 
  {
    for(int i=0; i<n; i++)
    {
      out[i] = R_NaN;
    }
  }
  
  return out;
}