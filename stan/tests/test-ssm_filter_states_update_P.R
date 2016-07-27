#function: ssm_filter_states_update_P
context("ssm_filter_states_update_P")
test_that("Stan function ssm_filter_states_update_P works", {
  f <- function(m, p, P, Z, Finv) {
    modfit <- test_stan_function("ssm_filter_states_update_P",
                                 data = list(m = m, p = p, P = P, Z = Z,
                                             Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  P <- rand_spd_mat(m)
  Z <- matrix(rnorm(m * p), p, m)
  Finv <- solve(rand_spd_mat(p))
  expected <- P - P %*% t(Z) %*% Finv %*% Z %*% P
  output <- f(m, p, P, Z, Finv)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})
