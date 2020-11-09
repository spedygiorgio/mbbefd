\name{mbbefd-package}
\alias{mbbefd-package}
\alias{mbbefd}
\docType{package}
\title{
\packageTitle{mbbefd}
}
\description{
The idea of this package emerged in 2013 from G.A. Spedicato who 
at this time worked in the area of quantitative risk assessment. 
In 2015, M. Gesmann and C. Dutang joined the project.
This project is hosted at \href{https://github.com/spedygiorgio/mbbefd}{github}.

This package contains the core functions of the two parametrizations
of the MBBEFD distribution
(distribution function, density, quantile functions, random generation, 
aka d, p, q, r)
as well as MBBEFD exposure curve (ec) and raw moments (m). 

This package also provides other distributions used for destruction rate
modelling, that is the beta, the shifted truncated Pareto 
and the generalized beta distributions.
Due to the presence of a total loss, a one-inflated version of the
previous distributions is also provided.

The vignette shows code snippets to fit the distribution to empirical data: \href{../doc/introduction_to_mbbefd.pdf}{Exposure rating, destruction rate models and the mbbefd package}.

}
\details{

\tabular{ll}{
Package: \tab mbbefd\cr
Type: \tab Package\cr
Version: \tab 0.8.9.1\cr
License: \tab GPL-2\cr
}


}
\author{
\packageAuthor{mbbefd}

Maintainer: \packageMaintainer{mbbefd}
}
\references{
Bernegger, Stefan. The Swiss Re Exposure Curves And The MBBEFD Distribution Class. ASTIN Bulletin (1997) 27:1, p99,
\doi{10.2143/AST.27.1.563208}.
}
\keyword{ package }
\seealso{
See \code{\link{mbbefd-distr}} for the MBBEFD distribution;\cr
see \code{\link{swissRe}}, \code{\link{exposureCurve}} for exposure curves;\cr
see \code{\link{gbeta}}, \code{\link{stpareto}} for finite-support distributions;\cr
see \code{\link{oidistribution}}, \code{\link{oibeta}}, \code{\link{oigbeta}}, \code{\link{oiunif}}, \code{\link{oistpareto}} for one-inflated distributions.
}
