#function: ssm_update_L_u
context("ssm_update_L_u")

test_ssm_update_L_u <- function(m, Z, K) {
  modfit <-
    test_stan_function("ssm_update_L_u",
                       data = list(m = m,
                                   Z = array(Z),
                                   K = array(K)))
  output <- rstan::extract(modfit, "output")[[1]]
  dim(output) <- dim(output)[-1]
  expected <- diag(m) - K %*% t(Z)
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  Z <- rand_vec(m)
  K <- as.numeric(rand_kalman_gain(m, 1))
  test_ssm_update_L_u(m, Z, K)
}
