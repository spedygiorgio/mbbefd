# `mbbefd` NEWS

## CHANGES IN `mbbefd` VERSION 0.8.16

(0.8.15 was a development version only; the last release on CRAN is 0.8.14, so the changes of both are listed here.)

### NEW

  * Allow parallel computation `bootDR()` as in `bootdist()`.  
  * Add a `control` argument for `fitDR()` which allows different levels of estimation traces.
  * Add components for `fitDR()` outputs in order to make correct bootstrap estimations in `bootDR()`.
  * Control arguments are correctly passed to optimization functions in `fitDR()`.
  * Add a `vcov()` function for `bootDR()` outputs.
  * Rewrite entirely the `d`, `p`, `q`, `r` functions for GB1 in order to allow a parameter vector.
  * Rewrite entirely `*oifun()` functions for one-inflated distributions in order to allow a vector `p1`.
  * Speedup `fitDR()` by removing `fitdist()` in favor of `mledist()` directly.
  * `fitDR()` for `mbbefd` and `MBBEFD` now explores the whole parameter space, including the region 
    `b < 1` and `gb > 1` (decreasing density, `a < -1` in the first parametrization), which was 
    previously excluded: the maximum likelihood is computed on three domains and the best fit is retained.
  * `dmbbefd()`, `pmbbefd()`, `qmbbefd()`, `rmbbefd()`, `ecmbbefd()`, `mmbbefd()` and `tlmbbefd()` now 
    support `a < -1` with `0 < b < 1`, consistently with the `MBBEFD(g,b)` functions (previously `NaN`).
  * `fitDR(method = "tlmme")` for `mbbefd` and `MBBEFD` now solves the total loss and mean equations exactly, 
    as in Bernegger (1997): `g = 1/etl(x)` and `b` is the unique root of the mean equation. 
    A warning is raised and `NA` is returned when `0 < etl(x) < mean(x) < 1` does not hold.
  * `fitDR(method = "mle")` for `mbbefd` and `MBBEFD` raises a warning when the data contain no total loss:
    the total loss probability is then estimated only from the shape of the distribution below 1.
  * `fitDR(method = "mle")` for `mbbefd` and `MBBEFD` refines each constrained fit with an unconstrained 
    maximization in `(log(g-1), log(b))`, where the log-likelihood is continuous across the domain boundaries 
    `b = 1` and `gb = 1`. When the MLE is close to `b = 1` (e.g. Swiss Re curve `c = 4`), the constrained 
    fits stopped short of the boundary (log-likelihood gaps up to 0.16 and `NaN` standard errors); 
    the maximum is now reached in both parametrizations.
  * `swissRe()` accepts a vector `c` and then returns a matrix with columns `b` and `g` (previously 
    a meaningless vector); the documentation now warns that `b` comes before `g` and that 
    parameters should be selected by name.
  * Add unit tests for gradients, constraint jacobians, parametrization consistency, `tlmme`, 
    standard errors and the new parameter domain.
  
### BUG FIX 

  * Fix an issue in `fitDR()` for distribution `oibeta` with total-loss-moment-matching estimation.
  * Fix the treatment of `NA` and `NaN` in probability distribution functions (`d`, `p`, `q`, `ec`, `m`, `tl`).
  * `bootDR()` now restores the user's `warn` and `show.error.messages` options even when the bootstrap loop fails or is interrupted.
  * The start-up message no longer prints an empty `Date:` line when the `Date` field is absent from `DESCRIPTION`.
  * Fix the gradient of the log-likelihood of `MBBEFD(g,b)` with respect to `b` (a factor `-2` was missing).
    Point estimates were mostly unaffected, but the Hessian, hence standard errors, correlation and 
    confidence intervals of `fitDR(x, "MBBEFD")`, were wrong: standard errors were too small 
    (coverage of 95% Wald intervals between 3% and 83% in simulations, now about 95%).
  * Fix the jacobians of the constraints `constrmbbefd1jac()` and `constrmbbefd2jac()`, and the constraint 
    `b > 0` in `constrmbbefd2()` (previously `a > 0` was checked twice); they were used by 
    `fitDR(x, "mbbefd", method = "tlmme")`.
  * `fitDR(x, "MBBEFD", method = "tlmme")` did not solve the moment equations (estimates could get stuck at `b = 1`).
  * `fitDR()` for `mbbefd` and `MBBEFD` now computes the Hessian (and standard errors) also when the 
    optimization converges in a single domain (previously `NA`).
  * Fix the starting value of `b` based on `F(1/2)` in `fitDR(x, "MBBEFD")`: the polynomial root is `sqrt(b)`, not `b`.
  * Fix the point mass at 1 returned by `dmbbefd()` and `dMBBEFD()` in the special case `gb = 1` (`a = Inf`): it is `b`, as given by `tlmbbefd()`, `pmbbefd()` and `rmbbefd()`, not `1 - b`.
  * Fix the standard errors of the shape parameters of one-inflated fits (`oibeta`, `oistpareto`, `oigbeta`) in `fitDR()`: they are estimated on the observations below 1 only, so their covariance matrix is now divided by the number of such observations instead of the total sample size (standard errors were too small by a factor `sqrt(1 - p1)`).
  * Fix `dgbeta(log = FALSE)` (hence `doigbeta()`) with vectorized parameters when some shape 
    parameters are invalid or some observations are outside `[0, 1]`: the data were indexed 
    with the wrong subset, giving wrong densities and recycling warnings.
  * Fix the exposure curves of one-inflated distributions (`ecoiunif()`, `ecoibeta()`, `ecoistpareto()`, 
    `ecoigbeta()`) when `x` contains values outside `[0, 1]` or missing values: the exposure curve of the 
    continuous part was indexed twice, so that the results were shifted (e.g. `ecoibeta(c(-0.1, 0), 2, 3, 0.1)` 
    returned 0.697 at 0 instead of 0).
  * `dmbbefd()` and `dMBBEFD()` now return the density at `x = 0` (it was 0, so that the log-likelihood of 
    data containing a zero destruction rate was `-Inf`).
  * `bootDR()` uses scalar logical operators (`&&`, `||`) in its argument checks.
  * Import `uniroot()` from `stats` (used by the exact `tlmme` estimator; `R CMD check` NOTE).
  * C++ random generators (`rmbbefd()`, `rMBBEFD()` with scalar parameters) now follow exactly the 
    parameter domains of the R generators: invalid or missing parameters give `NaN` (previously, 
    e.g., `rmbbefd(n, 0.5, NA)`, `rmbbefd(n, -1, 0.5)` or `rMBBEFD(n, 0.5, 0)` returned 1 and 
    `rmbbefd(n, -Inf, 0.5)` or `rmbbefd(n, NaN, 0.5)` returned random values), and degenerate 
    parameters no longer consume random numbers, as in R. For the same seed, the C++ and R generators 
    give the same variates and leave the same random number generator state.
  * `rmbbefd()` and `rMBBEFD()` raise the error `invalid arguments` for a negative, missing, infinite or 
    too large `n`, as base R (previously, the error came from the allocation inside the C++ code).
  * The R generators `rmbbefdR()` / `rMBBEFDR()` (used for vector parameters) return `NaN` for missing 
    parameters instead of failing, and `rmbbefdR()` no longer draws the case `a = Inf` twice.

### MISC

  * Raise the minimal required version of `fitdistrplus` to 1.2-2: `fitDR()` and `bootDR()` rely on the `calcvcov` argument and on `density.bootdist()` (introduced in 1.2-1; before: `unused argument (calcvcov = TRUE)`), and the standard errors of one-inflated fits assume that `mledist()` minimizes the average negative log-likelihood (since 1.2-2).
  * Remove the obsolete `src/Makevars*` files (`RCPP_USE_UNWIND_PROTECT` is the default of `Rcpp >= 0.12.18`) and the unneeded `SystemRequirements: GNU make`.
  * The `LICENCE` file and the README badge now state GPL-2, as `DESCRIPTION` (they said GPL-3).
  * C++ generators: constants (`1/log(b)`, thresholds) computed once outside the loops and the output 
    vector is not zero-filled before being written (`no_init`): 15-40% faster on the main cases. Checked 
    with valgrind and UndefinedBehaviorSanitizer (no error).
  * Code style checked with `lintr` (configuration in `.lintr`, new GitHub Actions workflow `lint`): 
    spaces around operators and after commas, `if (`, double quotes, `<-` for assignment, no semicolons, 
    lines of at most 120 characters, no trailing whitespace. These changes are whitespace and quoting 
    only: the parsed code is identical. The brace style of the package (opening brace on its own line) 
    and its object names (exported API) are kept, so `brace_linter`, `indentation_linter` and 
    `object_name_linter` are disabled, as well as `cyclocomp_linter` and `commented_code_linter`.
  * Revision of the vignette: corrected formulas (exposure curve of the beta and one-inflated beta 
    distributions, parameter domain of `MBBEFD(a,b)`, which includes `a < -1`), new sections on the 
    parameter domains, the Swiss Re curves, the likelihood with a point mass, the fitting algorithm 
    and the exact total-loss-moment matching. The premium rating example now follows Swiss Re (2004): 
    the risks above the retention of the surplus treaty are apportioned with the mean MPL (not the 
    maximum MPL), and the premium rate (1.65%, against 1.636% in the publication) is computed in the 
    vignette instead of being typed (the code returned `NaN`).
  * Comments added to the code and to the unit tests: each file starts with an overview of what it does, and 
    each function and test explains its role, its arguments and the formula it implements.

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

 
