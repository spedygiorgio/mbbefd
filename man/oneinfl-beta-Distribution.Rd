\name{oibeta}
\alias{oibeta}
\alias{doibeta}
\alias{poibeta}
\alias{qoibeta}
\alias{roibeta}
\alias{ecoibeta}
\alias{moibeta}

\title{
One-inflated beta distribution
}
\description{
These functions perform probabilistic analysis as well as random sampling 
on one-inflated beta distribution.

}
\usage{

doibeta(x, shape1, shape2, p1, ncp=0, log=FALSE)
poibeta(q, shape1, shape2, p1, ncp=0, lower.tail = TRUE, log.p = FALSE)
qoibeta(p, shape1, shape2, p1, ncp=0, lower.tail = TRUE, log.p = FALSE)
roibeta(n, shape1, shape2, p1, ncp=0)
ecoibeta(x, shape1, shape2, p1, ncp=0)
moibeta(order, shape1, shape2, p1, ncp=0)

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
\item{p1, shape1, shape2, ncp}{
parameters.
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
\code{d,p,q,ec,m}-\code{oibeta} functions computes the density function,
the distribution function, the quantile function, the exposure curve function and raw moments
of the one-inflated beta distribution.
\code{roibeta} generates random variates of this distribution.

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
#TODO
}

