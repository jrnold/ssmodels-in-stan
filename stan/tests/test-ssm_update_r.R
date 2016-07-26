test_that("Stan function ssm_smoth_update_r works", {
  f <- function(m, p, r, Z, v, Finv, L) {
    modfit <- test_stan_function("ssm_smooth_update_r",
                                 data = list(m = m, p = p, r = r, Z = Z,
                                             v = v, Finv = Finv, L = L))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  r <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  v <- rnorm(p)
  Finv <- solve(rand_spd_mat(p))
  L <- matrix(rnorm(m * m), m, m)
  expected <- t(Z) %*% Finv %*% v + t(L) %*% r
  output <- f(m, p, array(r), Z, array(v), Finv, L)
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})
