# `mbbefd` NEWS

## CHANGES IN `mbbefd` VERSION 0.8.15

### NEW

  * Allow parallel computation `bootDR()` as in `bootdist()`.  
  * Add a `control` argument for `fitDR()` which allows different levels of estimation traces.
  * Add components for `fitDR()` outputs in order to make correct bootstrap estimations in `bootDR()`.
  * Rewrite entirely the `d`, `p`, `q`, `r` functions for GB1 in order to allow a parameter vector.
  * Rewrite entirely `*oifun()` functions for one-inflated distributions in order to allow a vector `p1`.
  
### BUG FIX 

  * Fix an issue in `fitDR()` for distribution `oibeta` with total-loss-moment-matching estimation.

# `mbbefd` NEWS

## CHANGES IN `mbbefd` VERSION 0.8.14

  * Rewrite entirely the `d`, `p`, `q`, `r` functions for MBBEFD in order to allow a parameter vector.
  * Convert `NEWS` file in markdown.
  * Remove the use of `dilog()` from package `gsl`.
  * Remove a NOTE in R CMD CHECK with a missing link in a man page.

## CHANGES IN `mbbefd` VERSION 0.8.13

### NEW

 * Rewrite entirely the internal function `fitDR.addcomp()` to correctly compute asymptotic covariance matrix.
 * In `fitDR()`, the covariance matrix is now computed for MBBEFD (2nd parametrization) and is correctly computed for other distributions.
 * Remove some tests from the targz build.
 * Cleanup `NEWS` file.

## CHANGES IN `mbbefd` VERSION 0.8.12

### NOTE FIX

 * Remove missing links.

## CHANGES IN `mbbefd` VERSION 0.8.11

### NEW

 * Merge `ChangeLog` and `NEWS` files.
 * Now use Rd format for `NEWS` file.
 * Now use HTML2 for the vignette `test-beta`.
 * Add DOIs for some references.
 * Clean-up `NAMESPACE` file.

### TYPO

 * Remove package version in the overall man page.
 * Remove TODOs in man pages.

### WARNING FIX

 * Remove/update old URLs in man page.
 * Remove URLs in `README` using badger package.

## CHANGES IN `mbbefd` VERSION 0.8.10

### WARNING FIX

 * Fix a warning due to a missing link in documentation object from `fitDR.Rd`: to `fitdistrplus.quantile()`.

## CHANGES IN `mbbefd` VERSION 0.8.9.1

### NEW

 * Change the maintainer from G.A. Spedicato to C. Dutang.
 * Remove of the knitrcitations dependency.

### TYPO

 * Correct the vignette, in particular typo in Equation (10).

## CHANGES IN `mbbefd` VERSION 0.8-8.5

### NEW

 * Add ORCID.

### BUG FIX

 * For MBBEFD, the cdf now uses `res[q \<= 0] \<- 0` rather than `res[q \< 0] \<- 0`.

## CHANGES IN `mbbefd` VERSION 0.8-8.4

### NEW

 * Bumped requirement of `Rcpp` and `R`.
 * Added support for the unwind api of `Rcpp` [https://www.r-bloggers.com/2018/07/boost-the-speed-of-r-calls-from-rcpp/](https://www.r-bloggers.com/2018/07/boost-the-speed-of-r-calls-from-rcpp/).
 * Update automatic naming of legend in `eccomp()` (similar update for fitdistrplus).
 * Add a lines method for eecf object and update the man page.

### BUG FIX

 * Update testing parameter value of mbbefd(a,b) for inconsistent value.

## CHANGES IN `mbbefd` VERSION 0.8-8.3

### NEW

 * Update starting values to meet requirements of next fitdistrplus version.

### WARNING FIX

 * Remove unused arguments.

## CHANGES IN `mbbefd` VERSION 0.8-8.2

### WARNING FIX

 * Fix issues in docs.

## CHANGES IN `mbbefd` VERSION 0.8-8.2

### WARNING FIX

 * Removed vignette in inst... It will be autocreated.
 * Native routines registration fix.
 * Bumped requirements to R 3.4.

## CHANGES IN `mbbefd` VERSION 0.8-8

### NEW

 * Native routines registration (compliance with R 3.4).

## CHANGES IN `mbbefd` VERSION 0.8-7

### NEW

 * Minor changes (simplified integration).

## CHANGES IN `mbbefd` VERSION 0.8-6

### BUG FIX

 * fix in `fitDR()` and `doi(g)beta()`.
 * Add corresponding test on smokedfish dataset.

## CHANGES IN `mbbefd` VERSION 0.8-5

### WARNING FIX

 * Minor changes to cope with fitdistrplus changes.

## CHANGES IN `mbbefd` VERSION 0.8

### NEW

 * Added data sets (aon and loss).
 * Start an unified approach for paramater estimation based on the fitdist function from fitdistrplus.

## CHANGES IN `mbbefd` VERSION 0.7

### NEW

 * Markus Gesmann added as contributor, adding new vignette on exposure rating.

## CHANGES IN `mbbefd` VERSION 0.6.3

### BUG FIX

 * Fixed some calls to Rcpp code.

## CHANGES IN `mbbefd` VERSION 0.6.2

### NEW

 * Christophe Dutang added as contributor.
 * Introduced shifted truncated Pareto functions.

## CHANGES IN `mbbefd` VERSION 0.6.1

### NEW

 * Harmonize distribution parameters' to R usual ones: `pmbbefd()` accepts (q, a, b, g) and no longer (x, a, b, g).
 * Improved vignettes.

## CHANGES IN `mbbefd` VERSION 0.6.0/1

### NEW

 * Move to github [https://github.com/spedygiorgio/mbbefd](https://github.com/spedygiorgio/mbbefd).
 * Rcpp introduced.

## CHANGES IN `mbbefd` VERSION 0.5.0.1

### WARNING FIX

 * Inserted missing packages for LaTeX.

## CHANGES IN `mbbefd` VERSION 0.5.0.1

### BUG FIX

 * Fixed various bugs in code and vignettes.

## CHANGES IN `mbbefd` VERSION 0.1

### NEW

 * First release on CRAN.

 
