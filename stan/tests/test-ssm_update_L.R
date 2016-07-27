test_that("Stan function ssm_update_L works", {
  f <- function(m, p, Z, T, K) {
    modfit <- test_stan_function("ssm_update_L",
                                 data = list(m = m, p = p, Z = Z, T = T, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  m <- 3L
  Z <- matrix(rnorm(m * p), p, m)
  T <- matrix(rnorm(m * m), m, m)
  P <- rand_spd_mat(m)
  # Need Kalman gain to be generated in this way; otherwise (T - K Z) not symm
  K <- T %*% P %*% t(Z) %*% solve(rand_spd_mat(p))
  expected <- T - K %*% Z
  output <- f(m, p, Z, T, K)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})
