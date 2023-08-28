\name{stpareto}
\alias{stpareto}
\alias{dstpareto}
\alias{pstpareto}
\alias{qstpareto}
\alias{rstpareto}
\alias{mstpareto}
\alias{ecstpareto}

\title{
The shifted truncated Pareto distribution
}
\description{
These functions perform probabilistic analysis as well as random sampling on the shifted truncated Pareto distribution.

}
\usage{

dstpareto(x, a, log=FALSE)
pstpareto(q, a, lower.tail = TRUE, log.p = FALSE)
qstpareto(p, a, lower.tail = TRUE, log.p = FALSE)
rstpareto(n, a)
mstpareto(order, a)
ecstpareto(x, a)


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
\item{order}{
order of the raw moment.
}
\item{a}{
shape parameter. 
}
\item{log, log.p}{
logical; if \code{TRUE}, probabilities \code{p} are given as log(p).
}
\item{lower.tail}{
logical; if \code{TRUE} (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X> x]}.
}
  

}
\details{
The distribution is based on the Pareto 2 truncated at 1. The distribution function is given
by \eqn{P(X<=x) = (1-(x+1)^(-a))/(1-2^(-a)).}

}
\value{
A numeric value or a vector.
}

\author{
Dutang Christophe
}
\seealso{
\code{\link{mbbefd-distr}}, \code{\link{exposureCurve}}
}
\examples{

#density
curve(dstpareto(x, 3))

#cdf
curve(pstpareto(x, 3))
}
\keyword{distribution}
