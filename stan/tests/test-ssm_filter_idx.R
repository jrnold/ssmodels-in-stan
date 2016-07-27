#function: ssm_filter_idx
#function: ssm_filter_get_a
#function: ssm_filter_get_P
#function: ssm_filter_get_v
#function: ssm_filter_get_Finv
#function: ssm_filter_get_K
#function: ssm_filter_size
#function: ssm_filter_get_loglik

context("ssm_filter_idx")
test_that("Stan function ssm_filter_idx and ssm_filter_get_* functions work", {


  f <- function(m, p, x) {
    modfit <- test_stan_function("ssm_filter_get",
                                 data = list(m = m, p = p, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  symmat_size <- function(n) n * (n + 1) / 2
  xsz <- 1 + p + symmat_size(p) + m * p + m + symmat_size(m)
  x <- seq_len(xsz)
  output <- f(m, p, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(output[["idx"]][1, , ],
               matrix(c(1, 1, 1,
                      p, 1 + 1, 1 + p,
                      symmat_size(p), 2 + p, 1 + p + symmat_size(p),
                      m * p, 2 + p + symmat_size(p), 1 + p + symmat_size(p) + m * p,
                      m, 2 + p + symmat_size(p) + m * p, 1 + p + symmat_size(p) + m * p + m,
                      symmat_size(m), 2 + p + symmat_size(p) + m * p + m, 1 + p + symmat_size(p) + m * p + m + symmat_size(m)),
                      byrow = TRUE, 6, 3))
  expect_equivalent(output[["ll"]], array(1))
  expect_equivalent(output[["v"]][1, ], 2:5)
  expect_equivalent(output[["Finv"]][1, , ],
                    matrix(c(6, 7, 8, 9,
                             7, 10, 11, 12,
                             8, 11, 13, 14,
                             9, 12, 14, 15),
                           4, 4))
  expect_equivalent(output[["K"]][1, ,],
                    matrix(16:(16 + m * p - 1), m, p))
  expect_equivalent(as.numeric(output[["a"]]),
                    27 + 1:m)
  expect_equivalent(output[["P"]][1, ,],
                    matrix(c(31, 32, 33,
                             32, 34, 35,
                             33, 35, 36),
                           3, 3))
})
