test_that("Stan function ssm_filter_update_P works", {
  f <- function(m, p, P, Z, T, RQR, K) {
    modfit <- test_stan_function("ssm_filter_update_P",
                                 data = list(m = m, p = p, P = P,
                                             Z = Z, T = T, RQR = RQR, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  P <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  T <- matrix(rnorm(m * m), m, m)
  # Need Kalman gain to be generated in this way; otherwise (T - K Z) not symm
  K <- T %*% P %*% t(Z) %*% solve(rand_pdmat(p))
  RQR <- rand_pdmat(m)
  expected <- T %*% P %*% t(T - K %*% Z) + RQR
  output <- f(m, p, P, Z, T, RQR, K)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})
