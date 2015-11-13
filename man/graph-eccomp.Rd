\name{graph-eccomp}
\alias{graph-eccomp}
\alias{eccomp}

\title{Graphical comparison of multiple fitted distributions}
\description{   
 \code{eccomp} plots the empirical exposure curve distribution 
 against fitted exposure curve functions.
}

\usage{
eccomp(ft, xlim, ylim, main, xlab, ylab, do.points=TRUE,
                    datapch, datacol, fitlty, fitcol, addlegend = TRUE, 
                   legendtext, xlegend = "bottomright", 
                    ylegend = NULL, \dots)
    
}
\arguments{
\item{ft}{One \code{"DR"} object or a list of objects of class \code{"DR"}.}
\item{xlim}{The \eqn{x}-limits of the plot.}
\item{ylim}{The \eqn{y}-limits of the plot.}
\item{main}{A main title for the plot, see also \code{\link{title}}.}
\item{xlab}{A label for the \eqn{x}-axis, defaults to a description of \code{x}.}
\item{ylab}{A label for the \eqn{y}-axis, defaults to a description of \code{y}.}
\item{datapch}{An integer specifying a symbol to be used in plotting data points, 
    see also \code{\link{points}}.}
\item{datacol}{A specification of the color to be used in plotting data points.}
\item{fitcol}{A (vector of) color(s) to plot fitted distributions. 
    If there are fewer colors than fits they are recycled in the standard fashion.}
\item{fitlty}{A (vector of) line type(s) to plot fitted distributions/densities. 
    If there are fewer colors than fits they are recycled in the standard fashion.
    See also \code{\link{par}}.}
\item{addlegend}{If \code{TRUE}, a legend is added to the plot.}
\item{legendtext}{A character or expression vector of length \eqn{\geq 1} to appear
    in the legend, see also \code{\link{legend}}.}
\item{xlegend, ylegend}{The \eqn{x} and \eqn{y} co-ordinates to be used to position 
  the legend. They can be specified by keyword or in any way which is
          accepted by 'xy.coords': see \code{\link{legend}} for details.}
\item{do.points}{logical; if \code{TRUE}, also draw points at the x-locations. 
  Default is true. For large dataset (n > 1e4), \code{do.points} is ignored and no
  point is drawn.}
\item{\dots}{Further graphical arguments passed to graphical functions used in cdfcomp, denscomp,
    ppcomp and qqcomp.}
}

\details{
    \code{eccomp} provides a exposure curve plot of each fitted distribution
    along with the eecf. 
    
    
     
    By default a legend is added to these plots. Many graphical arguments are optional, 
    dedicated to personalize the plots, and fixed to default values if omitted.
}

\seealso{ 
    See \code{\link{plot}}, \code{\link{legend}}, \code{\link{eecf}}.  
}



\author{ 
Christophe Dutang.
}

\examples{
# (1) 

}

\keyword{ distribution }
