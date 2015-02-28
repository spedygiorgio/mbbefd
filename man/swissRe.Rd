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
BERNEGGER, STEFAN. THE SWISS RE EXPOSURE CURVES AND THE MBBEFD DISTRIBUTION CLASS. Astin Bulletin (1997): 99.
}
\author{
Giorgio Spedicato
}

\seealso{
\code{\link{dmbbefd}}
}
\examples{
pars<-swissRe(4)
losses<-rmbbefd(n=1000,b=pars[1],g=pars[2])
mean(losses)
}
