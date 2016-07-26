test_that("Stan function ssm_filter_update_Finv_miss works", {
  f <- function(m, p, p_t, Z, P, H, y_idx) {
    modfit <- test_stan_function("ssm_filter_update_Finv_miss",
                                 data = list(m = m, p = p, p_t = p_t,
                                   Z = Z, P = P,
                                   H = H, y_idx = y_idx))
    if (p_t == 0) {
      expected <- matrix(0, p, p)
    } else if (p_t == p) {
      expected <- solve(Z %*% P %*% t(Z) + H)
    } else {
      idx <- y_idx[1:p_t]
      Z_star <- Z[idx, , drop = FALSE]
      H_star <- H[idx, idx, drop = FALSE]
      expected <- matrix(0, p, p)
      Finv_star <- solve(Z_star %*% P %*% t(Z_star) + H_star)
      for (i in seq_along(idx)) {
        for (j in seq_along(idx)) {
          expected[idx[i], idx[j]] <- Finv_star[i, j]
        }
      }
    }
    output <- rstan::extract(modfit, "output")[[1]][1, , ]
    expect_length(output, p * p)
    expect_equal(output, expected, tolerance = 10e-5)
  }

  # No missing values
  m <- 3L
  p <- 4L
  p_t <- p
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_spd_mat(m)
  H <- rand_spd_mat(p)
  f(m, p, p_t, Z, P, H, y_idx = seq_len(p_t))

  # Some missing values
  p_t <- 2
  f(m, p, p_t, Z, P, H, y_idx = c(1, 3, 0, 0))

  # All missing values
  p_t <- 0
  f(m, p, p_t, Z, P, H, y_idx = rep(0, p))

})
