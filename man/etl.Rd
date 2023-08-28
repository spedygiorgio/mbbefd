\name{etl}
\alias{etl}

\title{
Empirical total loss
}
\description{
Compute the empirical total loss.

}
\usage{
etl(x, na.rm=FALSE)

}

\arguments{
\item{x}{numeric vector of the observations.}
\item{na.rm}{a logical value indicating whether \code{NA} values should 
be stripped before the computation proceeds.}

}
\details{
Compute the empirical total loss defined as the proportion of full destruction rates,
that is observations that equal 1.

}
\value{
A numeric value or a vector.
}
\author{
Dutang Christophe
}

\examples{
x <- c(1, 0.000495134903027804, 0.787229130724068, 0.71154311082138, 
0.0669802789251427, 0.310872967333683, 1, 1, 1, 1, 0.162030982251957, 
1, 1, 0.322530106394859, 1, 1, 1, 0.60805410798081, 0.660941675188664, 1)

#empirical total loss (true value is 1/2)
etl(x)
}

