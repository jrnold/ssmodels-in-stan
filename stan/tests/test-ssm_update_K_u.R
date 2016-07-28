#function: ssm_update_K_u
context("ssm_update_K_u")

test_ssm_update_K_u <- function(m, P, Z, Finv) {
  modfit <-
    test_stan_function("ssm_update_K_u",
                       data = list(m = m, P = P,
                                   Z = array(Z),
                                   T = T, Finv = Finv))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- as.numeric(P %*% Z %*% Finv)
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  P <- rand_spd_mat(m)
  Z <- rnorm(m)
  Finv <- exp(rnorm(1))
  test_ssm_update_K_u(m, P, Z, Finv)
}
