#function: ssm_smooth_eta_get_mean
#function: ssm_smooth_eta_get_mean
#function: ssm_smooth_eta_size

context("ssm_smooth_eta_get")
test_that("It works", {
  f <- function(q, x) {
    modfit <- test_stan_function("ssm_smooth_eta_get",
                                 data = list(q = q, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  q <- 3L
  xsz <- 9L
  x <- seq_len(xsz)
  output <- f(q, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["eta_mean"]]), 1:3)
  expect_equal(output[["eta_var"]][1, , ],
               matrix(c(4,  5,  6,
                        5,  7,  8,
                        6,  8,  9), 3, 3, byrow = TRUE))
})
