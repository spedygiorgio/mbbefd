\name{oiunif}
\alias{oiunif}
\alias{doiunif}
\alias{poiunif}
\alias{qoiunif}
\alias{roiunif}
\alias{ecoiunif}
\alias{moiunif}
\alias{tloiunif}

\title{
One-inflated uniform distribution
}
\description{
These functions perform probabilistic analysis as well as random sampling 
on one-inflated uniform distribution.

}
\usage{

doiunif(x, p1, log=FALSE)
poiunif(q, p1, lower.tail = TRUE, log.p = FALSE)
qoiunif(p, p1, lower.tail = TRUE, log.p = FALSE)
roiunif(n, p1)
ecoiunif(x, p1)
moiunif(order, p1)
tloiunif(p1)

}

\arguments{
\item{x, q}{
vector of quantiles.
}
\item{p}{
vector of probabilities.
}
\item{n}{
number of observations. If \code{length(n) > 1}, the length is take to be the number required.
}
\item{p1}{
parameter.
}
\item{order}{
order of the raw moment.
}
\item{log, log.p}{
logical; if \code{TRUE}, probabilities \code{p} are given as log(p).
}
\item{lower.tail}{
logical; if \code{TRUE} (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X> x]}.
}

}
\details{
\code{d,p,q,ec,m,tl}-\code{oiunif} functions computes the density function,
the distribution function, the quantile function, the exposure curve function, 
raw moments and total loss of the one-inflated uniform distribution.
\code{roiunif} generates random variates of this distribution.
}
\value{
A numeric value or a vector.
}

\author{
Dutang Christophe
}
\seealso{
\code{\link{mbbefd-distr}} and \code{\link{oidistribution}}.
}
\examples{

#density
curve(doiunif(x, 1/3), n=200, ylim=0:1)

#cdf
curve(poiunif(x, 1/3), n=200)
}
\keyword{distribution}
