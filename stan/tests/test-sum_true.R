#function: sum_true
context("sum_true")
test_that("Stan function sum_true works", {
  f <- function(x) {
    n <- length(x)
    modfit <- test_stan_function("int_sum_true",
                                 data = list(n = n, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(1, -1, 3, 0, 0)), 2L)
  expect_equal(f(c(1, 1, 1)), 3L)
  expect_equal(f(c(0, 0, 0)), 0L)

})
