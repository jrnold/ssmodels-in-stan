#function: ssm_update_a_u1
context("ssm_update_a_u1")

test_ssm_update_a_u1 <- function(m, a, K, v) {
  modfit <-
    test_stan_function("ssm_update_a_u1",
                       data = list(m = m,
                                   a = array(a),
                                   K = array(K),
                                   v = v))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- as.numeric(a + K %*% matrix(v))
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  y <- rnorm(1)
  a <- rand_vec(m)
  K <- rand_simplex(m)
  v <- rnorm(1)
  test_ssm_update_a_u1(m, a, K, v)
}
