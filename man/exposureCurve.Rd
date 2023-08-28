\name{exposureCurve}
\alias{exposureCurve}
\alias{ecbeta}
\alias{ecunif}

\title{
Exposure curves for the beta and the uniform distributions.
}
\description{
An exposure curve is defined between x between 0 and 1 and represents the ratio of the limited expected value to unlimited expected value.
}
\usage{
ecbeta(x, shape1, shape2)
ecunif(x, min = 0, max =1)

}
\arguments{

  \item{x}{
x value, percentage of damage to total loss
}
\item{shape1, shape2}{
parameters for the beta distribution.
}
\item{min, max}{
parameters for the uniform distribution.
}
}
\details{
\code{ecbeta}, \code{ecunif} is the theoretical exposure curve function for beta and uniform distribution.

}
\value{
A numeric value
}
\references{
BERNEGGER, STEFAN (1997). 
\emph{The Swiss Re Exposure Curves And The MBBEFD Distribution Class}, 
ASTIN Bulletin, 27(1), pp99-111, \doi{https://doi.org/10.2143/AST.27.1.563208}.
}
\author{
Giorgio Spedicato, 
Christophe Dutang
}
\seealso{
\code{ecmbbefd} and \code{ecMBBEFD} are implemented in \code{\link{mbbefd-distr}}. 
See also \code{\link{Uniform}}, \code{\link{Beta}}, \code{\link{swissRe}}.
}
\examples{

x <- 0.2
ecbeta(x, 2, 3)
ecunif(x)


}
