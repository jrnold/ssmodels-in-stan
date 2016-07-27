#function: stationary_cov
context("stationary_cov")

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

  # ARMA(1, 1)
  arma11 <- function(phi, theta) {
    x <- matrix(0, 2, 2)
    x[1, 1] <- (1 + theta ^ 2 + 2 * theta * phi) / (1 - phi ^  2)
    x[2, 1] <- x[1, 2] <- theta
    x[2, 2] <- theta ^ 2
    x
  }
  phi <- 0.5
  theta <- 0.25
  arma11(0.5, 0.25)
  T <- matrix(c(phi, 0, 1, 0), 2, 2)
  R <- matrix(c(1, theta))
  RQR <- tcrossprod(R)
  expect_equal(f(T, RQR), arma11(phi, theta))

  # AR(1, 1)
  phi <- 0.5
  T <- matrix(phi)
  RQR <- matrix(1)
  f(T, RQR)
  expect_equal(as.numeric(f(T, RQR)), 1 / (1 - phi ^ 2), tolerance = TOL)  
})
