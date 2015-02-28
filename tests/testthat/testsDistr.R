#unit testing
library(testthat)
library(mbbefd)

context("Check some values with Mahler")

test_that("some values", {
  expect_equal(round(mbbefdExposure(x=1/2,a = .2,b = .04),2), 0.68)
  expect_equal(1/dG(x=0,a=.2,b=.04),.6)
  expect_equal(round(dmbbefd(x = .1, a=.2, b=.04),2), .65)
  expect_equal(round(1-pmbbefd(q= .6, a=.2, b=.04),4),.5043)  
})

