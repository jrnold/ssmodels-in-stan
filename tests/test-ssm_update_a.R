#function: ssm_update_a
context("ssm_update_a")
test_that("Stan function ssm_update_a works", {
  f <- function(m, p, a, c, T, v, K) {
    modfit <- test_stan_function("ssm_update_a",
                                 data = list(m = m, p = p, a = a,
                                             c = c, T = T, v = v, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  a <- rnorm(m)
  c <- rnorm(m)
  T <- matrix(rnorm(m * m), m, m)
  v <- rnorm(p)
  K <- matrix(runif(m * p), m, p)
  expected <- T %*% a + K %*% v + c
  anew <- f(m, p, array(a), array(c), T, array(v), K)
  expect_length(anew, m)
  expect_equal(matrix(anew), expected, tolerance = 10e-5)
})
