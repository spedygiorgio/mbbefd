#unit tests for the fitting of the MBBEFD distribution (both parametrizations)
library(testthat)
library(mbbefd)

context("Fitting of the MBBEFD distribution")

#central finite differences (avoid a dependency on numDeriv)
numgrad <- function(f, x, h=1e-6)
  vapply(seq_along(x), function(i) { e <- replace(numeric(length(x)), i, h*max(1, abs(x[i])))
    (f(x+e) - f(x-e)) / (2*e[i]) }, numeric(1))
numjac <- function(f, x, h=1e-6) #rows = components of f, columns = components of x
  matrix(vapply(seq_along(x), function(i) { e <- replace(numeric(length(x)), i, h*max(1, abs(x[i])))
    (f(x+e) - f(x-e)) / (2*e[i]) }, numeric(length(f(x)))), ncol=length(x))
numhess <- function(f, x) numjac(function(t) numgrad(f, t, h=1e-5), x, h=1e-5)

test_that("analytical gradient of the log-likelihood", {
  set.seed(1)
  for(z in list(list("mbbefd", c(.5, .1)), list("mbbefd", c(-.5, 3)), list("mbbefd", c(-4/3, .5)),
                list("MBBEFD", c(5, 2)), list("MBBEFD", c(5, .1)), list("MBBEFD", c(5, .5))))
  {
    th <- z[[2]]
    x <- if(z[[1]] == "mbbefd") rmbbefd(200, th[1], th[2]) else rMBBEFD(200, th[1], th[2])
    ga <- mbbefd:::grLLfunc(x, th, z[[1]])
    gn <- numgrad(function(t) mbbefd:::LLfunc(x, t, z[[1]]), th)
    expect_equal(ga, gn, tolerance=1e-5)
  }
})

test_that("jacobians of the constraints", {
  for(z in list(list("constrmbbefd1", c(-.5, 3)), list("constrmbbefd2", c(.5, .1)),
                list("constrmbbefd3", c(-2, .5)), list("constrMBBEFD1", c(5, 2)),
                list("constrMBBEFD2", c(5, .1)), list("constrMBBEFD3", c(5, .5))))
  {
    hin <- get(z[[1]], envir=asNamespace("mbbefd"))
    jac <- get(paste0(z[[1]], "jac"), envir=asNamespace("mbbefd"))
    expect_equal(jac(z[[2]]), numjac(function(t) hin(t), z[[2]]), tolerance=1e-6)
  }
})

test_that("both parametrizations agree in the region b < 1, gb > 1 (a < -1)", {
  g <- 5; b <- 1/2; a <- (g-1)*b/(1-g*b)
  x <- c(0, .1, .5, .9, 1); p <- c(0, .2, .5, .79, .81, 1)
  expect_equal(dmbbefd(x, a, b), dMBBEFD(x, g, b))
  expect_equal(pmbbefd(x, a, b), pMBBEFD(x, g, b))
  expect_equal(qmbbefd(p, a, b), qMBBEFD(p, g, b))
  expect_equal(ecmbbefd(x, a, b), ecMBBEFD(x, g, b))
  expect_equal(mmbbefd(1:2, a, b), mMBBEFD(1:2, g, b))
  expect_equal(tlmbbefd(a, b), tlMBBEFD(g, b))
  set.seed(2)
  expect_false(any(is.nan(rmbbefd(100, a, b))))
  expect_false(any(is.nan(rmbbefd(100, rep(a, 2), rep(b, 2)))))
})

test_that("tlmme solves exactly the mean and total loss equations", {
  set.seed(3)
  for(th in list(c(5, .1), c(10, 3), c(5, .5)))
  {
    x <- rMBBEFD(500, th[1], th[2])
    f <- fitDR(x, "MBBEFD", method="tlmme")
    expect_equal(unname(1/f$estimate["g"]), etl(x), tolerance=1e-10)
    expect_equal(mMBBEFD(1, f$estimate["g"], f$estimate["b"]), mean(x), tolerance=1e-8)
    fa <- fitDR(x, "mbbefd", method="tlmme")
    expect_equal(tlmbbefd(fa$estimate["a"], fa$estimate["b"]), etl(x), tolerance=1e-8)
    expect_equal(mmbbefd(1, fa$estimate["a"], fa$estimate["b"]), mean(x), tolerance=1e-8)
  }
  expect_warning(f <- fitDR(c(.1, .2, .3), "MBBEFD", method="tlmme"))
  expect_true(all(is.na(f$estimate)))
})

test_that("mle reaches the region b < 1, gb > 1", {
  set.seed(4)
  x <- rMBBEFD(1000, 5, 1/2)
  f <- fitDR(x, "MBBEFD")
  expect_true(f$estimate["b"] < 1 && prod(f$estimate) > 1)
  expect_equal(unname(f$estimate), c(5, 1/2), tolerance=0.2)
  #the fit is at least as good as the true parameter value
  expect_true(f$loglik >= sum(dMBBEFD(x, 5, 1/2, log=TRUE)))
  fa <- fitDR(x, "mbbefd")
  expect_true(fa$estimate["a"] < -1)
  expect_equal(fa$loglik, f$loglik, tolerance=1e-4)
})

test_that("standard errors match the observed information", {
  set.seed(5)
  for(th in list(c(5, .1), c(10, 3)))
  {
    x <- rMBBEFD(1000, th[1], th[2])
    f <- fitDR(x, "MBBEFD")
    H <- -numhess(function(t) sum(dMBBEFD(x, t[1], t[2], log=TRUE)), unname(f$estimate))
    expect_equal(unname(f$sd), sqrt(diag(solve(H))), tolerance=1e-3)
    #b is a common parameter: its standard error does not depend on the parametrization
    fa <- fitDR(x, "mbbefd")
    expect_equal(unname(fa$sd["b"]), unname(f$sd["b"]), tolerance=1e-2)
  }
  #analytical Hessian of MBBEFD(a,b)
  x <- rmbbefd(300, -4/3, 1/2)
  expect_equal(mbbefd:::heLLfunc(x, c(-4/3, 1/2), "mbbefd"),
               numjac(function(t) mbbefd:::grLLfunc(x, t, "mbbefd"), c(-4/3, 1/2)), tolerance=1e-5)
})

test_that("a warning is raised without total losses", {
  set.seed(6)
  x <- rMBBEFD(200, 5, .1)
  expect_warning(fitDR(x[x < 1], "MBBEFD"), "no total loss")
})

test_that("standard errors of one-inflated fits match the observed information", {
  set.seed(7)
  x <- roibeta(1000, 2, 5, 1/10)
  f <- fitDR(x, "oibeta")
  y <- x[x < 1]
  H <- -numhess(function(t) sum(dbeta(y, t[1], t[2], log=TRUE)), unname(f$estimate[1:2]))
  expect_equal(unname(f$sd[1:2]), sqrt(diag(solve(H))), tolerance=1e-3)
  expect_equal(unname(f$sd[3]), sqrt(etl(x)*(1-etl(x))/length(x)), tolerance=1e-8)
})
