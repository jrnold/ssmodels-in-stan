#function: ssm_filter_states_get_a
#function: ssm_filter_states_get_P
#function: ssm_filter_states_size

context("ssm_filter_states_get")
test_that("Stan functions ssm_filter_states_get", {
  f <- function(m, x) {
    modfit <- test_stan_function("ssm_filter_states_get",
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
  expect_equal(as.numeric(output[["a"]]), 1:3)
  expect_equal(output[["P"]][1, , ],
               matrix(c(4, 5, 6,
                        5, 7, 8,
                        6, 8, 9), 3, 3))
})
