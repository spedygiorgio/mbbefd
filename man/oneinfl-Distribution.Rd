\name{oidistribution}
\alias{oidistribution}
\alias{doifun}
\alias{poifun}
\alias{qoifun}
\alias{roifun}
\alias{ecoifun}
\alias{moifun}

\title{
One-inflated distributions
}
\description{
These functions perform probabilistic analysis as well as random sampling 
on one-inflated distributions.

}
\usage{

doifun(x, dfun, p1, log=FALSE, ...)
poifun(q, pfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
qoifun(p, qfun, p1, lower.tail = TRUE, log.p = FALSE, ...)
roifun(n, rfun, p1, ...)
ecoifun(x, ecfun, mfun, p1, ...)
moifun(order, mfun, p1, ...)

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
\item{dfun, pfun, qfun, rfun}{
d, p, q, r functions of the original distribution. 
}
\item{p1}{
parameter for the probability at \code{x=1}.
}
\item{ecfun, mfun}{
exposure curve and moment functions which should have arguments
\code{x, ...} and \code{order, ...} respectively.
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
\item{\dots}{
further arguments to pass to \code{dfun, pfun, qfun, rfun, ecfun, mfun}.
}  

}
\details{
\code{d,p,q,ec,m} functions of \code{oifun} computes the density function,
the distribution function, the quantile function, the exposure curve function and raw moments
of an one-inflated distribution of an original distribution specified by \code{d,p,q,ec,m}-\code{fun}.
\code{roifun} generates random variates of the resulting distribution.

}
\value{
A numeric value or a vector.
}

\author{
Dutang Christophe
}
\seealso{
\code{\link{oibeta}}, \code{\link{oiunif}}, \code{\link{oistpareto}} and \code{\link{oidistribution}}.
}
\examples{
#TODO
}
\keyword{distribution}
