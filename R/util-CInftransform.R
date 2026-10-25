# File overview: smooth one-to-one maps between the real line and the intervals
# where the MBBEFD parameters live. They are used to turn a constrained
# maximum likelihood problem into an unconstrained one: the optimizer works on
# an unconstrained value u, and the parameter is obtained as Trans.*(u). The
# InvT.* functions give u back from a parameter value. They are used in
# fitDR-prefit.R to find starting values.

# Infinitely differentiable transformations of R into a bounded or half-bounded interval

#Transformation from (-Inf, +Inf) to (-1, 0)
# logistic map with a minus sign: u in R -> a in (-1, 0)
Trans.m10 <- function(x)
  -1/(1+exp(-x))
#Inverse
InvT.m10 <- function(x) 
  log(-x/(1+x))

#Transformation from (-Inf, +Inf) to (0, 1)
# logistic map: u in R -> b in (0, 1)
Trans.01 <- function(x)
  1/(1+exp(-x))
#Inverse
InvT.01 <- function(x)
  log(x/(1-x))

#Transformation from (-Inf, +Inf) to (0, +Inf)
# exponential map: u in R -> (0, +Inf)
Trans.0Inf <- function(x)
  exp(x)
#Inverse
InvT.0Inf <- function(x)
  log(x)

#Transformation from (-Inf, +Inf) to (1, +Inf)
# shifted exponential: u in R -> (1, +Inf), e.g. for g > 1 or b > 1
Trans.1Inf <- function(x)
  1+exp(x)
#Inverse
InvT.1Inf <- function(x)
  log(x-1)

#Transformation from (-Inf, +Inf) to (-Inf, -1)
# reflected shifted exponential: u in R -> (-Inf, -1), for a < -1 (region b < 1, gb > 1)
Trans.mInfm1 <- function(x)
  -1-exp(x)
#Inverse
InvT.mInfm1 <- function(x)
  log(-1-x)

