#function: ssm_stationary_cov
context("ssm_stationary_cov")
test_that("Stan function ssm_stationary_cov works", {
  f <- function(T, RQR) {
    modfit <- test_stan_function("stationary_cov",
                                 data = list(T = T, RQR = RQR,
                                             m = nrow(T)))
    ret <- rstan::extract(modfit)[["output"]]
    array(ret, dim(ret)[-1L])
  }

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
