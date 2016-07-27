#function: sum_false
context("sum_false")
test_that("Stan function sum_false works", {
  f <- function(x) {
    n <- length(x)
    modfit <- test_stan_function("int_sum_false",
                                 data = list(n = n, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(1, -1, 3, 0, 0)), 3L)
  expect_equal(f(c(1, 1, 1)), 0L)
  expect_equal(f(c(0, 0, 0)), 3L)

})
