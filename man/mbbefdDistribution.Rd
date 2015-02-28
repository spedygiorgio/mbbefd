\name{rmbbefd}
\alias{rmbbefd}
\alias{dmbbefd}
\alias{qmbbefd}
\alias{pmbbefd}
\title{
MBBEFD distribution definitions
}
\description{
These functions perform probabilistic analysis as well as random sampling on the MBBEFD distribution
}
\usage{

rmbbefd(n, a, b, g)

dmbbefd(x, a, b, g)

qmbbefd(p, a, b, g)

pmbbefd(q, a, b, g)
}

\arguments{

  \item{n}{
Sample size
}
  \item{a}{
a parameter
}
  \item{b}{
b parameter
}
  \item{g}{
g parameter (if a not given)
}

  \item{x}{
value of x
}

  \item{p}{
probability
}

  \item{q}{
quantile
}

}
\details{
it shall be remebered that \eqn{g=\frac{1}{p}=\frac{a+b}{\left(a+1\right)*b}}.
}
\value{
A numeric value or a vector.
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\author{
Giorgio Spedicato
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

