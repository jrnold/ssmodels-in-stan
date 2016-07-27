test_that("Stan function ssm_update_v_miss works", {
  f <- function(m, p, y, a, d, Z, p_t, y_idx) {
    modfit <- test_stan_function("ssm_update_v_miss",
                                 data = list(m = m, p = p, p_t = p_t,
                                   a = array(a),  d = array(d), y = array(y),
                                   Z = Z, y_idx = y_idx))
    expected <- rep(0, p)
    for (i in seq_len(p)) {
      if (i %in% y_idx) {
        expected[i] <- as.numeric(y[i] - d[i] - Z[i, ] %*% matrix(a))
      }
    }
    output <- as.numeric(rstan::extract(modfit, "output")[[1]][1, ])
    print(output)
    print(expected)
    expect_length(output, p)
    expect_equal(output, expected, tolerance = 10e-5)
  }

  # No missing values
  m <- 3L
  p <- 4L
  p_t <- p
  y <- rnorm(p)
  a <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  d <- rnorm(p)
  f(m, p, y, a, d, Z, p_t, y_idx = seq_len(p))

  # Some missing values
  p_t <- 2
  f(m, p, y, a, d, Z, p_t, y_idx = c(1, 3, 0, 0))

  # All missing values
  p_t <- 0
  f(m, p, y, a, d, Z, p_t, y_idx = rep(0, p))

})
