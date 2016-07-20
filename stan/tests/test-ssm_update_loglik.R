test_that("Stan function ssm_filter_update_ll works", {
  f <- function(m, p, v, Finv) {
    modfit <- test_stan_function("ssm_filter_update_ll",
                                 data = list(m = m, p = p, v = v, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    as.numeric(ret)
  }
  m <- 3L
  p <- 4L
  m <- 3L
  v <- rnorm(p)
  Finv <- solve(rand_pdmat(p))
  expected <- -0.5 * p * log(2 * base::pi) - 0.5 * (log(det(solve(Finv))) + t(v) %*% Finv %*% v)
  output <- f(m, p, array(v), Finv)
  expect_length(output, 1)
  expect_equal(output, as.numeric(expected), tol = 10e-5)
})
