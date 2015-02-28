\name{mbbefd-package}
\alias{mbbefd-package}
\alias{mbbefd}
\docType{package}
\title{
 MBBEFD distribution and exposure curve
}
\description{
It contains the definition of MBBEFD distribution (distribution function, density, quantile functions)
as well as MBBEFD exposure curve. The vignette shows code snippets to fit the distribution to empirical data.
}
\details{
\tabular{ll}{
Package: \tab mbbefd\cr
Type: \tab Package\cr
Version: \tab 0.6\cr
Date: \tab 2015-02-13\cr
License: \tab GPL 2\cr
}
}
\author{
Giorgio Spedicato
Maintainer:  Giorgio Spedicato <spedicato_giorgio@yahoo.it>
}
\references{
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\keyword{package}
\examples{
#sample variates
losses<-rmbbefd(n=10,a=0.2,b=0.04)
}
