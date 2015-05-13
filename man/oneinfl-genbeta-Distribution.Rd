\name{oigbeta}
\alias{oigbeta}
\alias{doigbeta}
\alias{poigbeta}
\alias{qoigbeta}
\alias{roigbeta}
\alias{ecoigbeta}
\alias{moigbeta}

\title{
One-inflated generalized beta of the first kind (GB1)) distribution
}
\description{
These functions perform probabilistic analysis as well as random sampling 
on one-inflated GB1 distribution.

}
\usage{

doigbeta(x, shape0, shape1, shape2, p1, log=FALSE)
poigbeta(q, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
qoigbeta(p, shape0, shape1, shape2, p1, lower.tail = TRUE, log.p = FALSE)
roigbeta(n, shape0, shape1, shape2, p1)
ecoigbeta(x, shape0, shape1, shape2, p1)
moigbeta(order, shape0, shape1, shape2, p1)

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
\item{p1, shape0, shape1, shape2}{
shape parameters.
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
\code{d,p,q,ec,m}-\code{oigbeta} functions computes the density function,
the distribution function, the quantile function, the exposure curve function and raw moments
of the one-inflated GB1 distribution.
\code{roigbeta} generates random variates of this distribution.

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
\keyword{distribution}
