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
fitDR(x, dist, method="mle", start=NULL, optim.method="default", \dots)
    
}

\arguments{
\item{x}{A numeric vector.}
\item{dist}{A character string \code{"name"} naming a distribution among
  \code{"oiunif"}, \code{"oistpareto"}, \code{"oibeta"}, \code{"oigbeta"}, 
  \code{"mbbefd"}, \code{"MBBEFD"}.}
\item{method}{A character string coding for the fitting method: 
    \code{"mle"} for 'maximum likelihood estimation', 
    \code{"tlmme"} for 'total-loss-moment matching estimation'.}
\item{start}{A named list giving the initial values of parameters 
    of the named distribution
    or a function of data computing initial values and returning a named list.
    This argument may be omitted (default) for some distributions for which reasonable 
    starting values are computed (see the 'details' section of  
    \code{\link[fitdistrplus]{mledist}}).}
\item{optim.method}{\code{"default"} or an optimization method to pass to \code{\link{optim}}.}

\item{\dots}{Further arguments to be passed to \code{"fitdist"}
    when \code{method != "tlmme"}.
    See \code{\link[fitdistrplus]{fitdist}} for details on parameter estimation.}
}

\details{

    The fitted distribution (\code{dist}) has its d, p, q, r functions defined in the
    man page: \code{\link{oiunif}}, \code{\link{oistpareto}}, \code{\link{oibeta}}, 
    \code{\link{oigbeta}}, \code{\link{mbbefd}}, \code{\link{MBBEFD}}.

    The two possible fitting methods are described below:
\describe{
    \item{When \code{method="mle"}}{
    Maximum likelihood estimation consists in maximizing the log-likelihood.
    A numerical optimization is carried out in \code{\link[fitdistrplus]{mledist}} via \code{optim}
    to find the best values (see \code{\link[fitdistrplus]{mledist}} for details).
    For one-inflated distributions, the probability parameter is estimated
    by a closed-form formula and other parameters use a two-optimization procedures.
    }
    \item{When \code{method="tlmme"}}{
    Total loss and moment matching estimation consists in equalizing theoretical and empirical 
    total loss as well as theoretical and empirical moments.
    The theoretical and the empirical moments are matched numerically, 
    by minimization of the sum of squared differences between observed and theoretical 
    quantities (see \code{\link[fitdistrplus]{mmedist}} for details).
    }
    
}
        
    For one-inflated distributions, 
    by default, direct optimization of the log-likelihood (or other criteria depending
    of the chosen method) is performed using \code{\link{optim}}, 
    with the "L-BFGS-B" method for distributions characterized by more than 
    one parameter and the "Brent" method for distributions characterized by only 
    one parameter. Note that when errors are raised by \code{optim}, it's a good 
    idea to start by adding traces during the optimization process by adding 
    \code{control=list(trace=1, REPORT=1)}.
    For the MBBEFD distribution, \code{\link[alabama]{constrOptim.nl}} is used.
    
    A pre-fitting process is carried out for the following distributions
    \code{"mbbefd"}, \code{"MBBEFD"} and \code{"oigbeta"} before
    the main optimization.
    
    The estimation process is carried out via \code{fitdist} from the 
    \code{fitdistrplus} package and the output object will inherit from the 
    \code{"fitdist"} class.
    Therefore, the following generic methods are available \code{print}, \code{plot}, 
    \code{summary}, \code{quantile}, \code{logLik}, \code{vcov} and \code{coef}.
    
    
}

\value{ 
    \code{fitDR} returns an object of class \code{"fitDR"} inheriting 
    from the \code{"fitdist"} class. That is a list with the following components:
    \item{estimate }{ the parameter estimates.}
    \item{method }{ the character string coding for the fitting method : 
        \code{"mle"} for 'maximum likelihood estimation', \code{"tlmme"} for 'matching total loss moment estimation'.}
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
        to functions \code{\link[fitdistrplus]{gofstat}}, \code{\link[fitdistrplus]{plotdist}} 
        and \code{\link[fitdistrplus]{cdfcomp}}. }
    \item{dots}{ the list of  further arguments passed in \dots to be used in \code{\link[fitdistrplus]{bootdist}} 
        in iterative calls to \code{\link[fitdistrplus]{mledist}}, \code{\link[fitdistrplus]{mmedist}},
        \code{\link[fitdistrplus]{qmedist}}, \code{\link[fitdistrplus]{mgedist}} or 
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
    \code{\link[fitdistrplus]{plotdist}}. An object of class "fitdist" or a list of objects of class
    "fitDR" corresponding to various fits using the same data set may also be plotted
    using a cdf plot (function \code{\link[fitdistrplus]{cdfcomp}}),
    a density plot(function \code{\link[fitdistrplus]{denscomp}}),
    a density Q-Q plot (function \code{\link[fitdistrplus]{qqcomp}}),
    or a P-P plot (function \code{\link[fitdistrplus]{ppcomp}}).
    }
    \item{\code{logLik}}{ Extracts the estimated log-likelihood from the \code{"fitDR"} object.
    }
    \item{\code{vcov}}{ Extracts the estimated var-covariance matrix from the 
      \code{"fitDR"} object (only available when \code{method = "mle"}).
    }
    \item{\code{coef}}{ Extracts the fitted coefficients from the \code{"fitDR"} object.
    }
}

}

\seealso{
    See \code{\link[fitdistrplus:mledist]{mledist}}, \code{\link[fitdistrplus:mmedist]{mmedist}}, 
    for details on parameter estimation.
    See \code{\link[fitdistrplus:gofstat]{gofstat}} for goodness-of-fit statistics.
    See \code{\link[fitdistrplus:plotdist]{plotdist}}, 
    \code{\link[fitdistrplus:graphcomp]{graphcomp}} for graphs.
    See \code{\link{bootDR}} for bootstrap procedures 
    See \code{\link{optim}} for base R optimization procedures.
    See \code{\link[fitdistrplus:quantile]{quantile.fitdist}}, another generic function, which calculates 
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

Delignette-Muller ML and Dutang C (2015), 
\emph{fitdistrplus: An R Package for Fitting Distributions}.
Journal of Statistical Software, 64(4), 1-34.

}


\author{ 
Christophe Dutang.
}

\examples{

# (1) fit of a one-inflated beta distribution by maximum likelihood estimation
#
n <- 1e3
set.seed(12345)
x <- roibeta(n, 3, 2, 1/6)

f1 <- fitDR(x, "oibeta", method="mle")
summary(f1)

plot(bootdist(f1, niter=11), enhance=TRUE, trueval=c(3, 2, 1/6))

}
\keyword{distribution}
