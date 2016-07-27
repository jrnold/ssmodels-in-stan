context("Stan function ssm_smooth_state_get")

test_that("they work", {
  f <- function(m, x) {
    modfit <- test_stan_function("ssm_smooth_state_get",
                                 data = list(m = m, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  xsz <- 9
  x <- seq_len(xsz)
  output <- f(m, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["alpha"]]), 1:3)
  expect_equal(output[["V"]][1, , ],
               matrix(c(4, 5, 6,
                        5, 7, 8,
                        6, 8, 9), 3, 3))
})
