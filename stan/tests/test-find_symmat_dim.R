context("find_symmat_dim")

test_that("Stan function find_symmat_dim works", {
  f <- function(n) {
    modfit <- test_stan_function("find_symmat_dim", data = list(n = n), init = NULL)
    as.numeric(rstan::extract(modfit, "output")[[1]])
  }
  expect_equal(f(1), 1)
  expect_equal(f(3), 2)
  expect_equal(f(6), 3)
})
