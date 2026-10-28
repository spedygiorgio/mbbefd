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
// The parameter domains, special cases and quantile formulas are the same as in
// the R generators rmbbefdR() and rMBBEFDR() (R/distr-mbbefdR-*.R): for the same
// seed, both give the same variates (up to rounding) and consume the same
// number of uniforms (none for degenerate or invalid parameters).
// Memory: the only buffer is the output vector, allocated by Rcpp and written
// at indices 0, ..., n-1; n is checked before the allocation, and errors are
// raised with Rcpp::stop() (a C++ exception converted into an R error by the
// code generated in RcppExports.cpp).
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// The number of variates must be a non-negative integer, as for the r*()
// functions of base R (NA_INTEGER is negative). Without this check, a negative
// length made Rf_allocVector() raise an R error from inside C++ code.
static void check_n(int n)
{
  if (n < 0)
    Rcpp::stop("invalid arguments");
}

// a vector of n identical values, without any random draw
static NumericVector constant_vector(int n, double value)
{
  NumericVector out = no_init(n);
  std::fill(out.begin(), out.end(), value);
  return out;
}


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
  
  check_n(n);
  // missing parameters: NaN
  if (ISNAN(a) || ISNAN(b))
    return constant_vector(n, R_NaN);
  // degenerate parameters: all the mass at 1 (total loss only), same cases as
  // in rmbbefdR(), no random draw
  if ((b == 0 && a > 0) || (a == -1 && b > 1) || (b == 1 && a != -1) ||
      (b == R_PosInf && a > -1 && a < 0) || (a == 0 && b > 0))
    return constant_vector(n, 1.0);
  
  NumericVector out = no_init(n); // every element is written below
  double u;
  // a = +Inf (gb = 1 in the 2nd param.): P(X = 1) = b and, below the total loss,
  // the quantile of the continuous part is log(1-U)/log(b)
  if (a == R_PosInf && b > 0 && b < 1)
  {
    // constants computed once (log() cannot be moved out of the loop by the
    // compiler, as it may set errno)
    const double thr = 1 - b, ilogb = 1 / std::log(b);
    for (int i = 0; i < n; i++)
    {
      u = unif_rand();
      out[i] = (u >= thr) ? 1.0 : std::log(1 - u) * ilogb;
    }
  // main domain: -1 < a < 0 with b > 1, a > 0 with 0 < b < 1, or a < -1 with
  // 0 < b < 1 (the last one corresponds to b < 1 and gb > 1 in the 2nd param.).
  // P(X = 1) = (a+1) b / (a + b) and the quantile of the continuous part is
  // log((1-U) a / (a + U)) / log(b).
  } else if (R_FINITE(a) && R_FINITE(b) &&
             ((-1 < a && a < 0 && b > 1) || (0 < a && b > 0 && b < 1) || (a < -1 && b > 0 && b < 1)))
  {
    const double thr = 1 - (a + 1) * b / (a + b), ilogb = 1 / std::log(b);
    for (int i = 0; i < n; i++)
    {
      u = unif_rand();
      out[i] = (u >= thr) ? 1.0 : std::log((1 - u) * a / (a + u)) * ilogb;
    }
  // invalid parameters
  } else
    return constant_vector(n, R_NaN);
  
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
  
  check_n(n);
  // missing parameters: NaN
  if (ISNAN(g) || ISNAN(b))
    return constant_vector(n, R_NaN);
  // degenerate parameters: all the mass at 1 (total loss only), same cases as
  // in rMBBEFDR(), no random draw
  if ((g == R_PosInf && b > 0 && b != 1) || (g > 1 && b == 0) || (g > 1 && b == R_PosInf) ||
      (g == 1 && b > 0 && b != 1))
    return constant_vector(n, 1.0);
  
  NumericVector out = no_init(n); // every element is written below
  double u;
  // special case gb = 1: P(X = 1) = b, X = log(1-U)/log(b) below the total loss
  if (g > 1 && b > 0 && b < 1 && b * g == 1)
  {
    const double thr = 1 - b, ilogb = 1 / std::log(b);
    for (int i = 0; i < n; i++)
    {
      u = unif_rand();
      out[i] = (u >= thr) ? 1.0 : std::log(1 - u) * ilogb;
    }
  // special case b = 1: X = U / ((1-U)(g-1)) below the total loss
  } else if (g > 1 && b == 1 && R_FINITE(g))
  {
    const double thr = 1 - 1 / g, gm1 = g - 1;
    for (int i = 0; i < n; i++)
    {
      u = unif_rand();
      out[i] = (u >= thr) ? 1.0 : u / (1 - u) / gm1;
    }
  // general case (b != 1, gb != 1): inverse of F(x) = 1 - (1-b)/((g-1) b^(1-x) + 1 - gb),
  // X = 1 - log((gb-1)/(g-1) + (1-b)/((1-U)(g-1))) / log(b) below the total loss
  } else if (g > 1 && b > 0 && b != 1 && b * g != 1 && R_FINITE(g) && R_FINITE(b))
  {
    const double thr = 1 - 1 / g, c1 = (g * b - 1) / (g - 1), c2 = (1 - b) / (g - 1);
    const double ilogb = 1 / std::log(b);
    for (int i = 0; i < n; i++)
    {
      u = unif_rand();
      out[i] = (u >= thr) ? 1.0 : 1 - std::log(c1 + c2 / (1 - u)) * ilogb;
    }
  // invalid parameters
  } else
    return constant_vector(n, R_NaN);
  
  return out;
}
