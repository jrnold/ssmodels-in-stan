#function: acf_to_pacf
context("acf_to_pacf")
test_that("Stan function acf_to_pacf works", {
  f <- function(x) {
    modfit <- test_stan_function("acf_to_pacf", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- ar_trans(rnorm(i))
    expect_equal(f(x), acf_to_pacf(x), tolerance = TOL)
  }
})
