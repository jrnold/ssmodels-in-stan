test_that("Stan function ssm_smoth_update_N works", {
  f <- function(m, p, N, Z, Finv, L) {
    modfit <- test_stan_function("ssm_smooth_update_N",
                                 data = list(m = m, p = p,
                                             N = N,
                                             Z = Z,
                                             Finv = Finv, L = L))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  N <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  Finv <- solve(rand_pdmat(p))
  L <- matrix(rnorm(m * m), m, m)
  expected <- t(Z) %*% Finv %*% Z + t(L) %*% N %*% L
  output <- f(m, p, N, Z, Finv, L)
  expect_equal(output, expected, tol = 10e-5)
})
