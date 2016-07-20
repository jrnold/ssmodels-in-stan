context("to_matrix_colwise")

test_that("Stan function to_matrix_colwise works", {
  f <- function(x, m, n) {
    modfit <- test_stan_function("to_matrix_colwise",
                       data = list(m = m, n = n, input = x))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  expect_equal(f(1:6, 2, 3), matrix(1:6, 2, 3))
})
