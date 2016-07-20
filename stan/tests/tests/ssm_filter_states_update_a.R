test_that("Stan function ssm_filter_states_update_a works", {
  f <- function(m, p, a, P = P, Z, v, Finv) {
    modfit <- test_stan_function("ssm_filter_states_update_a",
                                 data = list(m = m, p = p, a = a, P = P,
                                             Z = Z, v = v, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  a <- rnorm(m)
  P <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  v <- rnorm(p)
  Finv <- solve(rand_pdmat(p))
  expected <- a + P %*% t(Z) %*% Finv %*% v
  output <- f(m, p, a, P, Z, v, Finv)
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})
