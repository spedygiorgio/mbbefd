# File overview: the one-parameter family of Swiss Re exposure curves.
# Bernegger (1997, Section 4.3) showed that the Swiss Re curves Y1-Y4 and the
# Lloyd's curve correspond to MBBEFD(g,b) curves whose (g,b) lie on a smooth
# curve indexed by a single parameter c (c = 1.5, 2, 3, 4 for Y1-Y4, c = 5 for Lloyd's).

####Swiss Re curves####

# swissRe(c) returns the MBBEFD parameters c(b = ..., g = ...) of the Swiss Re
# curve with parameter c, using Bernegger's formulas
#   b(c) = exp(3.1 - 0.15 c (1 + c)),  g(c) = exp(c (0.78 + 0.12 c)).
# Note the order of the output: b first, then g, whereas the d,p,q,r,ec,...
# functions of the MBBEFD(g,b) parametrization take g before b: use the names
# (pars["g"], pars["b"]), not the positions. The order is kept for backward
# compatibility.
# For a scalar c the result is a named vector c(b, g); for a vector c it is a
# matrix with one row per value of c and columns "b" and "g" (before, a vector c
# returned a meaningless vector c(b1, b2, ..., g1, g2, ...) with wrong names).
swissRe <- function(c)
{
  stopifnot(is.numeric(c))
  b <- exp(3.1 - 0.15 * c * (1 + c))
  g <- exp(c * (0.78 + 0.12 * c))
  if (length(c) == 1)
    return(c(b = b, g = g))
  cbind(b = b, g = g)
}
