#empirical exposure curve
#in the spirit of ecdf
eecf <- function(x)
{
  f <- function(d)
    mean(pmin(x, d))/mean(x)
  Vectorize(f, "d")
}

#total loss
etl <- function(x)
  sum(x == 1)/length(x)

