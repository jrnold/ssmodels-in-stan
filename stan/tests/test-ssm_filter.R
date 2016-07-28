#function: ssm_filter
context("ssm_filter")

test_ssm_update_P_u2 <- function(m, P, .T, RQR) {
  modfit <-
    test_stan_function("ssm_update_P_u2",
                       data = list(m = m,
                                   P = P,
                                   T = .T,
                                   RQR = RQR))
  output <- rstan::extract(modfit, "output")[[1]]
  dim(output) <- dim(output)[-1]
  expected <- .T %*% P %*% t(.T) + RQR
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  P <- rand_spd_mat(m)
  .T <- rand_transition_mat(m)
  RQR <- rand_spd_mat(m)
  test_ssm_update_P_u2(m, P, .T, RQR)
}
