\name{mbbefdExposure}
\alias{mbbefdExposure}
\alias{dG}
\alias{ecbeta}
\alias{ecunif}

\title{
MBBEFD exposure curves and its derivative
}
\description{
The functions evaluate MBBEFD distribution exposure curve (and its derivative).
An exposure curve is defined between x between 0 and 1 and represents the ratio of the limited expected value to unlimited expected value
}
\usage{
mbbefdExposure(x, a, b, g)
dG(x, a, b, g)

ecbeta(x, shape1, shape2)
ecunif(x, min = 0, max =1)

}
\arguments{

  \item{x}{
x value, percentage of damage to total loss
}
  \item{a}{
a parameter
}
  \item{b}{
b parameter
}
  \item{g}{
g parameter, used if a not given.
}
\item{shape1, shape2}{
parameters for the beta distribution.
}
\item{min, max}{
parameters for the uniform distribution.
}
}
\details{
g parameter is the inverse of total loss probability and it is used if a not given.

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
}
\seealso{
\code{\link{dmbbefd}}, \code{\link{swissRe}}
}
\examples{
mbbefdExposure(x=0.2,b=0.04,g=20)
}
