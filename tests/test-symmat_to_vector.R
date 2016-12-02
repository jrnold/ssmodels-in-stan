#function: symmat_to_vector
context("symmat_to_vector")
test_that("Stan function symmat_to_vector works", {
  f <- function(x) {
    modfit <- test_stan_function("symmat_to_vector",
                              data = list(x = x,
                                          m = nrow(x),
                                          n = ncol(x)),
                              init = NULL)
    as.numeric(rstan::extract(modfit, "output")[[1]])
  }
  expect_equal(f(matrix(1, 1, 1)), 1)
  expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3)), as.numeric(1:6))
  # case with unequal rows/cols
  expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5), 3, 2)), c(1, 2, 4))
})
