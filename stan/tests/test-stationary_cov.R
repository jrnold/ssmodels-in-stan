context("Stan function stationary_cov")

test_that("Stan function ssm_loglik works", {
  f <- function(.T, .RQR) {
    .data <- list(m = nrow(.T),
                  T = .T,
                  RQR = .RQR)
    modfit <- test_stan_function("stationary_cov", data = .data)
    rstan::extract(modfit, "output")[[1]][1, , ]
  }
  m <- 3
  T <- rand_transition_mat(m)
  RQR <- rand_spd_mat(m)
  output <- f(T, RQR)
  expected <- matrix(solve(diag(m ^ 2) - kronecker(T, T), as.numeric(RQR)),
                     m, m)
  expect_equal(output, expected, TOL)
})
