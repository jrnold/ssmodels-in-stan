test_that("Stan function matrix_diff works", {
  f <- function(m, n, A, B) {
    modfit <- test_stan_function("matrix_diff",
                                 data = list(m = m, n = n,
                                             A = A, B = B))
    rstan::extract(modfit)[["output"]]
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  n <- 4L
  A <- matrix(rnorm(m * n), m, n)
  B <- matrix(rnorm(m * n), m, n)
  output <- f(m, n, A, B)
  expected <- max(abs(A - B)) / max(abs(A))
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})
