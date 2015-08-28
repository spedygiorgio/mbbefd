\name{fitDR}
\alias{fitDR}
\title{Fit of destruction rate models}

\description{
 Fit of univariate distributions to destruction rate data by maximum likelihood (mle), 
 moment matching (mme), quantile matching (qme) or 
 maximizing goodness-of-fit estimation (mge).
 The latter is also known as minimizing distance estimation.
 Generic methods are \code{print}, \code{plot}, 
\code{summary}, \code{quantile}, \code{logLik}, \code{vcov} and \code{coef}.
}

\usage{
fitDR(x, dist, method="mle", start=NULL, \dots)
    
}

\arguments{
\item{x}{A numeric vector.}
\item{dist}{A character string \code{"name"} naming a distribution among
  \code{"oiunif"}, \code{"oistpareto"}, \code{"oibeta"}, \code{"oigbeta"}, 
  \code{"mbbefd"}, \code{"MBBEFD"}.}
\item{method}{A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'moment matching estimation', 
    \code{"qme"} for 'quantile matching estimation' and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
\item{start}{A named list giving the initial values of parameters of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  
    \code{\link[fitdistrplus:fitdist]{mledist}}).}

\item{\dots}{Further arguments to be passed to \code{"fitdist"}.
    See \code{\link[fitdistrplus:fitdist]{fitdist}} for details on parameter estimation.}
}

\details{

    The fitted distribution (\code{dist}) has its d, p, q, r functions defined in the
    man page: \code{\link{oiunif}}, \code{\link{oistpareto}}, \code{\link{oibeta}}, \code{\link{oigbeta}}, 
    \code{\link{mbbefd}}, \code{\link{MBBEFD}}.

    The four possible fitting methods are described below:
\describe{
    \item{When \code{method="mle"}}{
    Maximum likelihood estimation consists in maximizing the log-likelihood.
    A numerical optimization is carried out in \code{\link{mledist}} via \code{optim}
    to find the best values (see \code{\link{mledist}} for details).
    }
    \item{When \code{method="mme"}}{
    Moment matching estimation consists in equalizing theoretical and empirical moments.
    The theoretical and the empirical moments are matched numerically, 
    by minimization of the sum of squared differences between observed and theoretical moments
    (see \code{\link{mmedist}} for details).
    }
    \item{When \code{method = "qme"}}{
    Quantile matching estimation consists in equalizing theoretical and empirical quantile.
    A numerical optimization is carried out in \code{\link{qmedist}} via \code{optim}
    to minimize of the sum of squared differences between observed and theoretical quantiles
    (see \code{\link{qmedist}} for details).
    }
    \item{When \code{method = "mge"}}{
    Maximum goodness-of-fit estimation consists in maximizing a goodness-of-fit statistics.
    A numerical optimization is carried out in \code{\link{mgedist}} via \code{optim}
    to minimize the goodness-of-fit distance (see \code{\link{mgedist}} for more details). 
    }
}
        
    By default, direct optimization of the log-likelihood (or other criteria depending
    of the chosen method) is performed using \code{\link{optim}}, 
    with the "Nelder-Mead" method for distributions characterized by more than one parameter
    and the "BFGS" method for distributions characterized by only one parameter. 
    Note that when errors are raised by \code{optim}, it's a good idea to start by adding traces during
    the optimization process by adding \code{control=list(trace=1, REPORT=1)}.
    
    The estimation process is carried out via \code{fitdist} from the \code{fitdistrplus} package
    and the output object will inherit from the \code{"fitdist"} class.
    Therefore, the following generic methods are available \code{print}, \code{plot}, 
    \code{summary}, \code{quantile}, \code{logLik}, \code{vcov} and \code{coef}.
    
    
}

\value{ 
    \code{fitDR} returns an object of class \code{"fitDR"} inheriting from the \code{"fitdist"} class.
    That is a list with the following components:
    \item{estimate }{ the parameter estimates.}
    \item{method }{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood estimation', \code{"mme"} for 'matching moment estimation',
         \code{"qme"} for 'matching quantile estimation'
    and \code{"mge"} for 'maximum goodness-of-fit estimation'.}
    \item{sd}{ the estimated standard errors, \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{cor}{ the estimated correlation matrix, \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{vcov}{ the estimated variance-covariance matrix, \code{NULL} if not available.}
    \item{loglik}{ the log-likelihood.}
    \item{aic}{ the Akaike information criterion.}
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion).}
    \item{n}{ the length of the data set.}
    \item{data}{ the data set.}
    \item{distname}{ the name of the distribution.}
    \item{fix.arg}{ the named list giving the values of parameters of the named distribution
        that must be kept fixed rather than estimated by maximum likelihood or \code{NULL} 
        if there are no such parameters. }
    \item{fix.arg.fun}{the function used to set the value of \code{fix.arg} or \code{NULL}.}    
    \item{discrete}{ the input argument or the automatic definition by the function to be passed
        to functions \code{\link{gofstat}}, \code{\link{plotdist}} 
        and \code{\link{cdfcomp}}. }
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{\link{bootdist}} 
        in iterative calls to \code{\link{mledist}}, \code{\link{mmedist}},
        \code{\link{qmedist}}, \code{\link{mgedist}} or 
        \code{NULL} if no such arguments.}
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}    

    Generic functions:
\describe{       
    \item{\code{print}}{
    The print of a \code{"fitDR"} object shows few traces about the fitting method and 
    the fitted distribution.
    }
    \item{\code{summary}}{
    The summary provides the parameter estimates of the fitted distribution, the log-likelihood, 
    AIC and BIC statistics and when the maximum likelihood is used, the standard errors of the 
    parameter estimates and the correlation matrix between parameter estimates.
    }
    \item{\code{plot}}{
    The plot of an object of class "fitDR" returned by \code{fitdist} uses the function 
    \code{\link{plotdist}}. An object of class "fitdist" or a list of objects of class
    "fitDR" corresponding to various fits using the same data set may also be plotted
    using a cdf plot (function \code{\link{cdfcomp}}),
    a density plot(function \code{\link{denscomp}}),
    a density Q-Q plot (function \code{\link{qqcomp}}),
    or a P-P plot (function \code{\link{ppcomp}}).
    }
    \item{\code{logLik}}{ Extracts the estimated log-likelihood from the \code{"fitDR"} object.
    }
    \item{\code{vcov}}{ Extracts the estimated var-covariance matrix from the \code{"fitDR"} object
    (only available When \code{method = "mle"}).
    }
    \item{\code{coef}}{ Extracts the fitted coefficients from the \code{"fitDR"} object.
    }
}

}

\seealso{
    See \code{\link[fitdistrplus:mledist]{mledist}}, \code{\link[fitdistrplus:mmedist]{mmedist}}, 
    \code{\link[fitdistrplus:qmedist]{qmedist}}, \code{\link[fitdistrplus:mgedist]{mgedist}} 
    for details on parameter estimation.
    See \code{\link[fitdistrplus:gofstat]{gofstat}} for goodness-of-fit statistics.
    See \code{\link[fitdistrplus:plotdist]{plotdist}}, 
    \code{\link[fitdistrplus:graphcomp]{graphcomp}} for graphs.
    See \code{\link{bootDR}} for bootstrap procedures 
    See \code{\link{optim}} for base R optimization procedures.
    See \code{\link[fitdistrplus:quantiles]{quantile.fitdist}}, another generic function, which calculates 
    quantiles from the fitted distribution.
    See \code{\link{quantile}} for base R quantile computation.
}

\references{ 
Cullen AC and Frey HC (1999), \emph{Probabilistic techniques in exposure assessment}.
Plenum Press, USA, pp. 81-155.

Venables WN and Ripley BD (2002), \emph{Modern applied statistics with S}.
Springer, New York, pp. 435-446.

Vose D (2000), \emph{Risk analysis, a quantitative guide}.
John Wiley & Sons Ltd, Chischester, England, pp. 99-143.

Delignette-Muller ML and Dutang C (2015), \emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Christophe Dutang.
}

\examples{

# (1) fit of a gamma distribution by maximum likelihood estimation
#


}
\keyword{distribution}
