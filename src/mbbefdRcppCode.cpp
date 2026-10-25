// File overview: C++ code compiled with Rcpp.
// - g2a(): conversion of the MBBEFD parameters (g,b) into a.
// - rmbbefdC() / rMBBEFDC(): random generation of MBBEFD(a,b) and MBBEFD(g,b)
//   by inversion: draw U uniform on (0,1); if U falls in the top probability
//   P(X = 1) the variate is a total loss (1), otherwise it is the quantile
//   function of the continuous part evaluated at U.
// The functions marked [[Rcpp::export]] are made callable from R by
// Rcpp::compileAttributes(), which (re)generates src/RcppExports.cpp and
// R/RcppExports.R: do not edit those two files by hand. The random generators
// are exported with a leading dot (.rmbbefdC, .rMBBEFDC) and called from
// R/distr-mbbefdCpp.R. unif_rand() uses R's random number generator, so
// set.seed() makes the draws reproducible.
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
// g2a: a = (g - 1) b / (1 - g b), the first parametrization from the second one.
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

// rmbbefdC: n variates of MBBEFD(a,b).
// [[Rcpp::export(.rmbbefdC)]]
NumericVector rmbbefdC(int n, double a, double b) {
  
  NumericVector out(n);
  double u, pab;
  // degenerate parameters: all the mass is at 1 (total loss only)
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
  // a = +Inf: P(X = 1) = 1 - b and X = log(1-U)/log(b) otherwise
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
  // main domain: -1 < a < 0 with b > 1, a > 0 with 0 < b < 1, or a < -1 with
  // 0 < b < 1 (the last one corresponds to b < 1 and gb > 1 in the 2nd param.).
  // P(X = 1) = (a+1) b / (a + b) and the quantile of the continuous part is
  // log((1-U) a / (a + U)) / log(b).
  }else if( (a +1 >0 && b > 0 && a*(1-b) > 0) || (a < -1 && b > 0 && b < 1) ) //main domain incl. b < 1 and gb > 1
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
  // invalid parameters
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

// rMBBEFDC: n variates of MBBEFD(g,b), with P(X = 1) = 1/g.
// [[Rcpp::export(.rMBBEFDC)]]
NumericVector rMBBEFDC(int n, double g, double b) {
  
  NumericVector out(n);
  double u;
  
  // degenerate parameters: total loss only
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
  // special case gb = 1: X = log(1-U)/log(b) below the total loss
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
  // special case b = 1: X = U / ((1-U)(g-1)) below the total loss
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
  // general case (b != 1, gb != 1): inverse of F(x) = 1 - (1-b)/((g-1) b^(1-x) + 1 - gb)
  }else if(b > 0 && g > 1)
  {
    for(int i=0; i<n; i++)
    {
      u = unif_rand();
      if(u > 1-1/g)
        out[i] = 1.0;
      else
        out[i] = 1-log((g*b-1)/(g-1) + (1-b)/((1-u)*(g-1)))/log(b);
    }
  // invalid parameters
  }else 
  {
    for(int i=0; i<n; i++)
    {
      out[i] = R_NaN;
    }
  }
  
  return out;
}