test_that("Stan functions ssm_filter_eta_get works", {
  f <- function(q, x) {
    modfit <- test_stan_function("ssm_smooth_eta_get",
                                 data = list(q = q, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  q <- 3L
  p <- 4L
  xsz <- 9
  x <- seq_len(xsz)
  output <- f(q, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["eta_mean"]]), 1:3)
  expect_equal(output[["eta_var"]][1, , ],
               matrix(c(4, 5,  6,
                        5, 7,  8,
                        6, 8, 9), 3, 3, byrow = TRUE))
})
