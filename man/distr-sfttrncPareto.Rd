\name{rstpareto}
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
TODO!

}
\value{
A numeric value or a vector.
}
\references{
TODO
}
\author{
Dutang Christophe
}
\seealso{
\code{\link{mbbefd-distr}}, \code{\link{exposureCurve}}
}
\examples{

dstpareto(0:4/4, 2)

pstpareto(0:4/4, 1/2)

}
\keyword{distribution}
