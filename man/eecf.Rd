\name{eecf}
\alias{eecf}
\alias{plot.eecf}
\alias{lines.eecf}
\alias{print.eecf}
\alias{summary.eecf}

\title{
Empirical Exposure Curve Function
}
\description{
Compute an empirical exposure curve function, with
several methods for plotting, printing, computing with
such an object.

}
\usage{

eecf(x)


\method{plot}{eecf}(x, \dots, ylab="Gn(x)", do.points=TRUE, 
     col.01line = "gray70", pch = 19, main=NULL, ylim=NULL, 
     add=FALSE)

\method{lines}{eecf}(x, \dots)

\method{print}{eecf}(x, digits= getOption("digits") - 2, \dots)

\method{summary}{eecf}(object, \dots)
}

\arguments{
  \item{x, object}{numeric vector of the observations for \code{eecf}; 
  for the methods, an object of class \code{"eecf"}.
  }
  \item{\dots}{arguments to be passed to subsequent methods, e.g.,
    to the \code{plot} method.}
  \item{ylab}{label for the y-axis.}
  \item{do.points}{logical; if \code{TRUE}, also draw points at 
  the (\code{xlim} restricted) knot locations.}
  \item{col.01line}{numeric or character specifying the color of the
    horizontal lines at y = 0 and 1, see \code{\link{colors}}.}
  \item{pch}{plotting character.}
  \item{main}{main title.}
  \item{ylim}{the y limits of the plot.}
  \item{add}{logical; if \code{TRUE} add to an already existing plot.}  
  \item{digits}{number of significant digits to use, see
    \code{\link{print}}.}

}
\details{
Compute a continuous empirical exposure curve and returns an object
of class \code{"eecf"} similar to what an object returned 
by \code{\link{ecdf}}.

}
\value{
For \code{eecf}, a function of class \code{"eecf"}, inheriting 
from the \code{"function"} class.

For the \code{summary} method, a summary of the knots of object 
with a \code{"header"} attribute.
}

\author{
Dutang Christophe
}
\seealso{
\code{\link{exposureCurve}}, \code{\link{ecdf}}.
}
\examples{

x <- c(0.4756816, 0.1594636, 0.1913558, 0.2387725, 0.1135414, 0.7775612,
  0.6858736, 0.4340655, 0.3181558, 0.1134244)

#print
eecf(x)

#summary
summary(eecf(x))

#plot
plot(eecf(x))

#lines
lines(eecf(x[1:4]), col="red")
}

