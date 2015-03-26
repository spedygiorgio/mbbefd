\name{mbbefd-distr}
\alias{mbbefd-distr}
\alias{rmbbefd}
\alias{dmbbefd}
\alias{qmbbefd}
\alias{pmbbefd}

\alias{dmbbefdR}
\alias{pmbbefdR}
\alias{qmbbefdR}
\alias{rmbbefdR}
\alias{ecmbbefdR}
\alias{mmbbefdR}
\alias{tlmbbefdR}

\alias{dMBBEFDR}
\alias{pMBBEFDR}
\alias{qMBBEFDR}
\alias{rMBBEFDR}
\alias{ecMBBEFDR}
\alias{mMBBEFDR}
\alias{tlMBBEFDR}

\title{
The MBBEFD distribution (two parametrizations)
}
\description{
These functions perform probabilistic analysis as well as random sampling on the MBBEFD distribution.

}
\usage{

rmbbefd(n, a, b, g)

dmbbefd(x, a, b, g)

qmbbefd(p, a, b, g)

pmbbefd(q, a, b, g)


dmbbefdR(x, a, b, log=FALSE)
pmbbefdR(q, a, b, lower.tail = TRUE, log.p = FALSE)
qmbbefdR(p, a, b, lower.tail = TRUE, log.p = FALSE)
rmbbefdR(n, a, b)
ecmbbefdR(x, a, b)
mmbbefdR(order, a, b)
tlmbbefdR(a, b)

dMBBEFDR(x, g, b, log=FALSE)
pMBBEFDR(q, g, b, lower.tail = TRUE, log.p = FALSE)
qMBBEFDR(p, g, b, lower.tail = TRUE, log.p = FALSE)
rMBBEFDR(n, g, b)
ecMBBEFDR(x, g, b)
mMBBEFDR(order, g, b)
tlMBBEFDR(g, b)
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
\item{a, b, g}{
shape parameters. For \code{.mbbefd} functions, \code{g} is computed from \code{a}.
}
\item{log, log.p}{
logical; if \code{TRUE}, probabilities \code{p} are given as log(p).
}
\item{lower.tail}{
logical; if \code{TRUE} (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X> x]}.
}
  

}
\details{
it shall be remebered that \eqn{g=\frac{1}{p_1}=\frac{a+b}{\left(a+1\right)*b}}.
}
\value{
A numeric value or a vector.
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\author{
Giorgio Spedicato,
Dutang Christophe
}
\seealso{
\code{\link{swissRe}}, \code{\link{mbbefdExposure}}
}
\examples{
#few examples of the functions
aPar=0.2
bPar=0.04
rmbbefd(n=1000,a=aPar,b=bPar) #for random generation
qmbbefd(p=0.7,a=aPar,b=bPar) #for quantiles
dmbbefd(x=0.5,b=.04,g=20) #for density
pmbbefd(q=0.5,b=.04,g=20) #for distribution function
}

