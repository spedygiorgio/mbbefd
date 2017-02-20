\name{mbbefd-package}
\alias{mbbefd-package}
\alias{mbbefd}
\docType{package}
\title{
 MBBEFD distribution, destruction rate and exposure curve modelling
}
\description{

The idea of this package emerged in 2013 from G Spedicato who at this time worked in the
area of quantitative risk assessment. In 2015, M Gesmann and C Dutang joined the project.
This project is hosted in \href{https://github.com/spedygiorgio/mbbefd}{github}.

This package contains the definition of the MBBEFD distributions
(distribution function, density, quantile functions, random generation, aka d, p, q, r)
as well as MBBEFD exposure curve (ec) and raw moments (m). 

This package also focuses on other distributions used for destruction rate
modelling, that is the beta, the shifted truncated Pareto 
and the generalized beta distributions.
Due to the presence of total loss, a one-inflated version of the
previous distributions is also provided.

The vignette shows code snippets to fit the distribution to empirical data: \href{../doc/introduction_to_mbbefd.pdf}{Exposure rating, destruction rate models and the mbbefd package}.
}
\details{
\tabular{ll}{
Package: \tab mbbefd\cr
Type: \tab Package\cr
Version: \tab 0.8-8\cr
Date: \tab 2017-02-22\cr
License: \tab GPL 2\cr
}
}
\author{
Giorgio Alfredo Spedicato [cre, aut], Christophe Dutang [aut]
and Markus Gesmann [ctb] \cr
Maintainer:  Giorgio Spedicato <spedicato_giorgio@yahoo.it>
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\keyword{package}
\seealso{
See \code{\link{mbbefd-distr}} for the MBBEFD distribution;\cr
see \code{\link{swissRe}}, \code{\link{exposureCurve}} for exposure curves;\cr
see \code{\link{gbeta}}, \code{\link{stpareto}} for finite-support distributions;\cr
see \code{\link{oidistribution}}, \code{\link{oibeta}}, \code{\link{oigbeta}}, \code{\link{oiunif}}, \code{\link{oistpareto}} for one-inflated distributions.
}

