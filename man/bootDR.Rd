\name{bootDR}
\alias{bootDR}

\title{ Bootstrap simulation of destruction rate models}

\description{
 Uses parametric or nonparametric bootstrap resampling in order to simulate uncertainty 
 in the parameters of the distribution fitted to destruction rate data.  
 Generic methods are \code{print}, \code{plot}, \code{summary}.
}

\usage{

bootDR(f, bootmethod="param", niter=1001, silent=TRUE)

}

\arguments{
\item{f}{ An object of class \code{"fitDR"}, output of the \code{\link{fitDR}} function.}
\item{bootmethod}{ A character string coding for the type of resampling : \code{"param"} 
  for a parametric resampling and \code{"nonparam"} for a nonparametric resampling of data.}
\item{niter}{ The number of samples drawn by bootstrap.  }
\item{silent}{A logical to remove or show warnings and errors when bootstraping.}

}

\details{
    Samples are drawn by parametric bootstrap (resampling from the distribution fitted by 
    \code{\link{fitDR}}) or nonparametric bootstrap (resampling with replacement from the 
    data set). On each bootstrap sample the estimaiton process is 
    used to estimate bootstrapped values of parameters. When that function fails
    to converge, \code{NA} values are returned. Medians and 2.5 and 97.5 percentiles are 
    computed by removing \code{NA} values. 
    
    This method returns an object of class \code{"bootDR"} inheriting from the 
    \code{"bootdist"} class. Therefore the following generic methods are defined:
     \code{print}, \code{plot}, \code{summary}.
    
    
}

\value{ 
    \code{bootDR} returns an object of class \code{"bootDR"} inheriting from the \code{"bootdist"} class.
    That is a list with 6 components,
    \item{estim}{ a data frame containing the bootstrapped values of parameters.}
    \item{converg}{ a vector containing the codes for convergence obtained if an iterative method is
        used to estimate parameters on 
        each bootstraped data set (and 0 if a closed formula is used).}
    \item{method}{ A character string coding for the type of resampling : \code{"param"} for a parametric
        resampling and \code{"nonparam"} for a nonparametric resampling. }
    \item{nbboot}{ The number of samples drawn by bootstrap.}   
    \item{CI}{ bootstrap medians and 95 percent confidence percentile intervals of parameters. }
    \item{fitpart}{ The object of class \code{"fitDR"} on which the bootstrap procedure
        was applied.}
        
Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"bootDR"} object shows the bootstrap parameter estimates. 
    If inferior to the whole number of bootstrap iterations, 
    the number of iterations for which the estimation converges is also printed.
    }
    \item{\code{summary}}{
    The summary provides the median and 2.5 and 97.5 percentiles of each parameter. 
    If inferior to the whole number of bootstrap iterations, 
    the number of iterations for which the estimation converges is also printed in the summary.
    }
    \item{\code{plot}}{
    The plot shows the bootstrap estimates with \code{\link{stripchart}} function
    for univariate parameters and \code{\link{plot}} function for multivariate parameters.
    }           
    
}
        
}
        

\seealso{ 
    See \code{\link[fitdistrplus:mledist]{mledist}}, \code{\link[fitdistrplus:mmedist]{mmedist}}, 
    \code{\link[fitdistrplus:qmedist]{qmedist}}, \code{\link[fitdistrplus:mgedist]{mgedist}} 
    for details on parameter estimation.
    See \code{\link[fitdistrplus:bootdist]{bootdist}} for details on generic function.
    See \code{\link{fitDR}} for estimation procedures.
    

}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 181-241.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.
}


\author{ 
Christophe Dutang
}


\examples{
# We choose a low number of bootstrap replicates in order to satisfy CRAN running times
# constraint.
# For practical applications, we recommend to use at least niter=501 or niter=1001.



}


\keyword{distribution}