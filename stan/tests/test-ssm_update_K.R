#function: ssm_update_K
context("ssm_update_K")
test_that("Stan function ssm_update_K works", {
  f <- function(m, p, P, Z, T, Finv) {
    modfit <- test_stan_function("ssm_update_K",
                                 data = list(m = m, p = p, T = T, Z = Z,
                                             P = P, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  T <- matrix(rnorm(m * m), m, m)
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_spd_mat(m)
  Finv <- solve(rand_spd_mat(p))
  expected <- T %*% P %*% t(Z) %*% Finv
  output <- f(m, p, P, Z, T, Finv)
  expect_length(output, m * p)
  expect_equal(output, expected, tolerance = 10e-5)
})
