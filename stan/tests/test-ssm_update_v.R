test_that("Stan function ssm_update_v works", {
  f <- function(m, p, y, a, d, Z) {
    modfit <- test_stan_function("ssm_update_v",
                                 data = list(m = m, p = p,
                                   a = array(a), d = array(d), y = array(y),
                                   Z = Z))

    expected <- as.numeric(y - d - Z %*% matrix(a))
    output <- as.numeric(rstan::extract(modfit, "output")[[1]][1, ])
    expect_length(output, p)
    expect_equal(output, expected, tolerance = 10e-5)
  }
  p <- 4L
  m <- 3L
  y <- rnorm(p)
  a <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  d <- rnorm(p)
  f(m, p, y, a, d, Z)
})
