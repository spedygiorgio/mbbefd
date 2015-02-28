\name{mbbefdExposure}
\alias{mbbefdExposure}
\alias{dG}

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
}
\details{
g parameter is the inverse of total loss probability and it is used if a not given.
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
