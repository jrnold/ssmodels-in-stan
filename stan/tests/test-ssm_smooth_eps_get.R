#function: ssm_smooth_eps_get
context("ssm_smooth_eps_get")
test_that("Stan functions ssm_filter_eps_get works", {
  f <- function(p, x) {
    modfit <- test_stan_function("ssm_smooth_eps_get",
                                 data = list(p = p, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  xsz <- 14
  x <- seq_len(xsz)
  output <- f(p, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["eps_mean"]]), 1:4)
  expect_equal(output[["eps_var"]][1, , ],
               matrix(c(5, 6,  7,  8,
                        6, 9,  10, 11,
                        7, 10, 12, 13,
                        8, 11, 13, 14), 4, 4, byrow = TRUE))
})
