#function: ssm_update_P_u1
context("ssm_update_P_u1")

test_ssm_update_P_u1 <- function(m, P, Finv, K) {
  modfit <-
    test_stan_function("ssm_update_P_u1",
                       data = list(m = m,
                                   P = P,
                                   Finv = Finv,
                                   K = array(K)))
  output <- rstan::extract(modfit, "output")[[1]]
  dim(output) <- dim(output)[-1]
  expected <- (P - K  %*% t(K) / Finv)
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  P <- rand_spd_mat(m)
  Finv <- exp(rnorm(1))
  K <- as.numeric(rand_kalman_gain(m, 1))
  test_ssm_update_P_u1(m, P, Finv, K)
}
