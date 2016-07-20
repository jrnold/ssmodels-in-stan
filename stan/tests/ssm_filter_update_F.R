test_that("Stan function ssm_filter_update_F works", {
  f <- function(m, p, Z, P, H) {
    modfit <- test_stan_function("ssm_filter_update_F",
                                 data = list(m = m, p = p, Z = Z, P = P, H = H))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_pdmat(m)
  H <- rand_pdmat(p)
  expected <- Z %*% P %*% t(Z) + H
  output <- f(m, p, Z, P, H)
  expect_length(output, p * p)
  expect_equal(output, expected, tolerance = 10e-5)
})
