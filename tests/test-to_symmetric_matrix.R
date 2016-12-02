#function: to_symmetric_matrix
context("to_symmetric_matrix")

test_that("Stan function to_symmetric_matrix works", {
  f <- function(x) {
    modfit <- test_stan_function("to_symmetric_matrix",
                       data = list(m = nrow(x), n = ncol(x), input = x))
    y <- rstan::extract(modfit, "output")[[1]]
    array(y, dim(y)[-1L])
  }
  expect_equal(f(matrix(1)), matrix(1))
  expect_equal(f(matrix(c(1, 2, 3, 4), 2, 2)),
               matrix(c(1, 2.5, 2.5, 4), 2, 2))
})
