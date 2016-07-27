#function: ssm_sim_idx
#function: ssm_sim_get_y
#function: ssm_sim_get_a
#function: ssm_sim_get_eta
#function: ssm_sim_get_eps
#function: ssm_sim_size

context("ssm_sim_get")
test_that("Stan functions ssm_sim_get works", {
  f <- function(m, p, q, x) {
    modfit <- test_stan_function("ssm_sim_get",
                                 data = list(m = m, p = p, q = q, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 5L;
  p <- 4L;
  q <- 3L;
  xsz <- 2 * p + m + q
  x <- seq_len(xsz)
  output <- f(m, p, q, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(output[["idx"]][1, , ],
               matrix(c(4, 1, 4,
                        5, 5, 9,
                        4, 10, 13,
                        3, 14, 16), 4, 3, byrow = TRUE))
  expect_equal(as.numeric(output[["y"]]), 1:4)
  expect_equal(as.numeric(output[["a"]]), 5:9)
  expect_equal(as.numeric(output[["eps"]]), 10:13)
  expect_equal(as.numeric(output[["eta"]]), 14:16)

})
