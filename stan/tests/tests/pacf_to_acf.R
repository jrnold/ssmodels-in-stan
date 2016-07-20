test_that("Stan function pacf_to_acf works", {
  f <- function(x) {
    modfit <- test_stan_function("pacf_to_acf", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- tanh(rnorm(i))
    expected <- pacf_to_acf(x)
    output <- f(x)
    expect_equal(output, expected, tolerance = TOL)
  }
})
