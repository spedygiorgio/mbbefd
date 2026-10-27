# File overview: exposure curve functions for the uniform and beta distributions.
# The exposure curve of a loss X on [0,1] is G(x) = E[min(X, x)] / E[X]: the share
# of the expected loss retained with a deductible x (as a fraction of the sum insured).
# The limited expected value E[min(X, x)] (lev*) and the moments (m*) come from actuar.
#new exposure curve function for some distributions
#limited expected value is implemented in actuar

# Exposure curve of the uniform distribution on [min, max].
ecunif <- function(x, min = 0, max = 1)
{
  levunif(x, min = min, max = max) / munif(1, min = min, max = max)
}

# Exposure curve of the beta distribution.
ecbeta <- function(x, shape1, shape2)
{
  levbeta(x, shape1 = shape1, shape2 = shape2) / mbeta(1, shape1 = shape1, shape2 = shape2)
}
