test_that("Stan function cholesky_decompose2 works", {
  f <- function(A) {
    n <- nrow(A)
    modfit <- test_stan_function("cholesky_decompose2",
                                 data = list(n = n, A = A))
    rstan::extract(modfit)[["output"]][1, , ]
  }
  A <- rand_pdmat(2)
  cholA <- t(chol(A))
  expect_equal(f(A), cholA, tolerance = 1e-5)
  A1 <- A[1, 1]
  cholA1 <- cholA[1, 1]
  expect_equal(f(matrix(c(A1, 0, 0, 0), 2, 2)),
                 matrix(c(cholA[1, 1], 0, 0, 0), 2, 2), tolerance = 1e-5)
  expect_equal(f(matrix(0, 2, 2)), matrix(0, 2, 2), tolerance = 1e-5)
})
