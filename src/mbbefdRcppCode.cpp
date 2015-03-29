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