#function: ssm_update_v_u
context("ssm_update_v_u")

test_ssm_update_v_u <- function(m, y, a, d, Z) {
  modfit <-
    test_stan_function("ssm_update_v_u",
                       data = list(m = m,
                                   y = y,
                                   a = array(a),
                                   Z = array(Z),
                                   d = d))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- as.numeric(y - Z %*% a - d)
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  y <- rnorm(1)
  a <- rnorm(m)
  d <- rnorm(1)
  Z <- rnorm(m)
  test_ssm_update_v_u(m, y, a, d, Z)
}
