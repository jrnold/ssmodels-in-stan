test_that("Stan function kronecker_prod works", {
  f <- function(A, B) {
    modfit <- test_stan_function("kronecker_prod",
                              data = list(A = A, B = B,
                                          m = nrow(A), n = ncol(A),
                                          p = nrow(B), q = ncol(B)),
                              init = NULL)
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  A <- matrix(1:2, 1, 2)
  B <- matrix(3 + 1:12, 3, 4)
  expect_equal(f(A, B), kronecker(A, B))
})
