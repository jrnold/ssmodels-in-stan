#function: unconstrain_stationary
context("unconstrain_stationary")
test_that("Stan function unconstrain_stationary works", {
  f <- function(x) {
    modfit <- test_stan_function("unconstrain_stationary", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    # real numbers
    expected <- rnorm(i)
    # AR coefficients
    phi <- ar_trans(expected)
    output <- f(phi)
    expect_equal(f(phi), expected, tolerance = TOL)
  }
})
