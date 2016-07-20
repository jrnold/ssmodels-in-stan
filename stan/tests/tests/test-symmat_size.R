context("symmat_size")

test_that("Stan function symmat_size works", {
  f <- function(n) {
    modfit <- test_stan_function("symmat_size", data = list(n = n), init = NULL)
    ret <- rstan::extract(modfit, "output")[[1]]
    as.numeric(ret)
  }
  expect_equal(f(3), 6)
  expect_equal(f(2), 3)
  expect_equal(f(1), 1)
})
