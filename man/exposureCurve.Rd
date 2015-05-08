\name{exposureCurve}
\alias{exposureCurve}
\alias{ecbeta}
\alias{ecunif}

\title{
Exposure curves for the MBBEFD, the beta and the uniform distributions.
}
\description{
An exposure curve is defined between x between 0 and 1 and represents the ratio of the limited expected value to unlimited expected value
}
\usage{
ecmbbefd(x, a, b)
ecMBBEFD(x, g, b)

ecbeta(x, shape1, shape2)
ecunif(x, min = 0, max =1)

}
\arguments{

  \item{x}{
x value, percentage of damage to total loss
}
\item{a, b, g}{
shape parameters for the MBBEFD distribution. 
}
\item{shape1, shape2}{
parameters for the beta distribution.
}
\item{min, max}{
parameters for the uniform distribution.
}
}
\details{
\code{ecmbbefd} and \code{ecMBBEFD} corresponds to the two parametrization of MBBEFD, see \code{\link{mbbefd-distr}} for details.

\code{ecbeta}, \code{ecunif} is the theoretical exposure curve function for beta and uniform distribution.

}
\value{
A numeric value
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\author{
Giorgio Spedicato
Christophe Dutang
}
\seealso{
\code{\link{mbbefd-distr}}, \code{\link{Uniform}}, \code{\link{Beta}}, \code{\link{swissRe}}.
}
\examples{

x <- 0.2
ecmbbefd(x, 2, 0.04)
ecMBBEFD(x, 20, 0.04)

ecbeta(x, 2, 3)
ecunif(x)


}
