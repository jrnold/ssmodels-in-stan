#function: constrain_stationary
context("constrain_stationary")
test_that("Stan function constrain_stationary works", {
  f <- function(x) {
    modfit <- test_stan_function("constrain_stationary", data = list(n = length(x), x = array(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- rnorm(i)
    output <- f(x)
    expected <- ar_trans(x)
    expect_equal(output, expected, tolerance = 10e-5)
  }
})
