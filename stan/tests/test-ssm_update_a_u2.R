#function: ssm_update_a_u2
context("ssm_update_a_u2")

test_ssm_update_a_u2 <- function(m, a, .c, .T) {
  modfit <-
    test_stan_function("ssm_update_P_u1",
                       data = list(m = m,
                                   a = array(a),
                                   c = array(.c),
                                   T = .T))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- as.numeric(.c + .T %*% a)
  expect_equal(output, expected, TOL)
}

for (m in c(1, 2, 3)) {
  a <- rand_vec(m)
  .c <- rand_vec(m)
  .T <- rand_transition_mat(m)
  test_ssm_update_a_u2(m, a, .c, .T)
}
