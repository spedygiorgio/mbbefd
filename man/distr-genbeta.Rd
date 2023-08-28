\name{gbeta}
\alias{gbeta}
\alias{dgbeta}
\alias{pgbeta}
\alias{qgbeta}
\alias{rgbeta}
\alias{ecgbeta}
\alias{mgbeta}

\alias{dgbeta1}


\title{The generalized Beta of the first kind Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the GB1 distribution with parameters \code{shape0}, 
  \code{shape1} and \code{shape2}.
}
\usage{
dgbeta(x, shape0, shape1, shape2, log = FALSE)
pgbeta(q, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
qgbeta(p, shape0, shape1, shape2, lower.tail = TRUE, log.p = FALSE)
rgbeta(n, shape0, shape1, shape2)
ecgbeta(x, shape0, shape1, shape2)
mgbeta(order, shape0, shape1, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape0, shape1, shape2}{positive parameters of the GB1 distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
    \item{order}{
order of the raw moment.
}

}
\details{
  The GB1 distribution with parameters \code{shape0} \eqn{= g}, 
  \code{shape1} \eqn{= a} and  \code{shape2} \eqn{= b} has density
  \deqn{f(x)=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}{x}^{a/g-1} {(1-x^{1/g})}^{b-1}/g%
  }{\Gamma(a+b)/(\Gamma(a)\Gamma(b))x^(a/g-1)(1-x^{1/g})^(b-1)/g}
  for \eqn{a,b,g > 0} and \eqn{0 \le x \le 1}
  where the boundary values at \eqn{x=0} or \eqn{x=1} are defined as
  by continuity (as limits).
  
}
\value{
  \code{dgbeta} gives the density, \code{pgbeta} the distribution
  function, \code{qgbeta} the quantile function, and \code{rgbeta}
  generates random deviates.

}

\references{
  Becker, R. A., Chambers, J. M. and Wilks, A. R. (1988)
  \emph{The New S Language},  Wadsworth & Brooks/Cole,
  \doi{10.1201/9781351074988}.

  Abramowitz, M. and Stegun, I. A. (1972)
  \emph{Handbook of Mathematical Functions.} New York: Dover.
  Chapter 6: Gamma and Related Functions.

  Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995)
  \emph{Continuous Univariate Distributions}, Volume 2, especially
  Chapter 25. Wiley, New York, \doi{10.1080/00224065.1996.11979675}.
}
\seealso{
  \link{Distributions} for other standard distributions.
}
\examples{

#density
curve(dgbeta(x, 3, 2, 3))

#cdf
curve(pgbeta(x, 3, 2, 3))
}
\keyword{distribution}
