\name{oistpareto}
\alias{oistpareto}
\alias{doistpareto}
\alias{poistpareto}
\alias{qoistpareto}
\alias{roistpareto}
\alias{ecoistpareto}
\alias{moistpareto}

\title{
One-inflated shifted truncated pareto distribution
}
\description{
These functions perform probabilistic analysis as well as random sampling 
on one-inflated shifted truncated pareto distribution.

}
\usage{

doistpareto(x, a, p1, log=FALSE)
poistpareto(q, a, p1, lower.tail = TRUE, log.p = FALSE)
qoistpareto(p, a, p1, lower.tail = TRUE, log.p = FALSE)
roistpareto(n, a, p1)
ecoistpareto(x, a, p1)
moistpareto(order, a, p1)

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
\item{a, p1}{
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
\code{d,p,q,ec,m}-\code{oistpareto} functions computes the density function,
the distribution function, the quantile function, the exposure curve function and raw moments
of the one-inflated shifted truncated pareto distribution.
\code{roistpareto} generates random variates of this distribution.

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

