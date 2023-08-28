\name{swissRe}
\alias{swissRe}

\title{
Swiss Re exposure curve generation function
}
\description{
This function turns out the MBBEFD b and g parameters for the famous Swiss Re (SR) exposure curves.
}
\usage{
swissRe(c)
}

\arguments{
  \item{c}{A numeric value}
}
\details{
The four Swiss Re Y1-Y4 are defined for c=1.5, 2, 3, 4. In addition c=5 coincides with a curve used by Lloyds for industrial risks exposure rating. 
}
\value{
A named two dimensional vector
}
\references{
BERNEGGER, STEFAN (1997). 
\emph{The Swiss Re Exposure Curves And The MBBEFD Distribution Class}, 
ASTIN Bulletin, 27(1), pp99-111, \doi{https://doi.org/10.2143/AST.27.1.563208}.
}
\author{
Giorgio Spedicato
}

\seealso{
\code{\link{mbbefd-distr}}.
}
\examples{
pars <- swissRe(4)
losses <- rMBBEFD(n=1000,b=pars[1],g=pars[2])
mean(losses)
}
