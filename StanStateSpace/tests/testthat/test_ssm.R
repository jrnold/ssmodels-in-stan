context("test ssm functions")

test_that("ssm_extract throws errors for invalid p, q, m arguments", {
  p <- 2
  m <- 4
  q <- 3
  iter <- 3
  time <- 2
  veclen <- 1 + p + p * (p + 1) / 2 + m * p + m + m * (m + 1) / 2
  x <- aperm(array(rep(seq_len(veclen * time), iter),
             c(veclen, time, iter)))
  expect_error(ssm_extract(x, m = 0, p = p, q = q),
               regexp = "m is not a count")
  expect_error(ssm_extract(x, m = c(1, 3), p = p, q = q),
               regexp = "m is not a count")
  expect_error(ssm_extract(x, m = m, p = 0, q = q),
               regexp = "p is not a count")
  expect_error(ssm_extract(x, m = m, p = 1:2, q = q),
               regexp = "p is not a count")
  expect_error(ssm_extract(x, m = m, p = p, q = 0),
               regexp = "q is not a count")
  expect_error(ssm_extract(x, m = m, p = p, q = 1:2),
               regexp = "q is not a count")
  expect_error(ssm_extract(x, m = m, p = p, q = 5),
               regexp = "q not less than or equal to m")
  expect_error(ssm_extract(array(x, dim = prod(dim(x))),
                           m = m, p = p, q = q),
               regexp = "not equal to 3L")
  expect_error(ssm_extract(as.character(x), m = m, p = p, q = q),
               regexp = "x is not a numeric or integer vector")
  expect_error(ssm_extract(x, m = m, p = p, q = q, type = "zzz"),
               regexp = "'arg' should be one of")
  expect_error(ssm_extract(x[ , , 1:2], m = m, p = p, q = q),
               regexp = "expected dim\\(x\\)\\[3\\]")
  expect_error(ssm_extract(x, m = m, p = p, q = q,
                           params = c("zzz")),
               regexp = "param.*value.*are invalid")

})

test_that("ssm_extract works for filter", {
  m <- 4
  p <- 2
  # v, Finv, K, a, P
  veclen <- 1 + p + p * (p + 1) / 2 + m * p + m + m * (m + 1) / 2
  time <- 5
  iter <- 6
  x <- array(rep(seq_len(veclen), each = time * iter), c(iter, time, veclen))
  out <- ssm_extract(x, m, p, type = "filter")
  expect_true(all(out$loglik ==  1))
  expect_equal(out$v[1, 1, ], c(2, 3))
  expect_equal(out$Finv[1, 1, , ], matrix(c(4, 5, 5, 6), 2, 2))
  expect_equal(out$K[1, 1, , ], matrix(6 + seq_len(m * p), m, p))
  expect_equal(out$a[1, 1, ], as.numeric(15:18))
  expect_equal(out$P[1, 1, , ], matrix(c(19, 20, 21, 22,
                                         20, 23, 24, 25,
                                         21, 24, 26, 27,
                                         22, 25, 27, 28), 4, 4, byrow = TRUE))
})

test_that("ssm_extract works for sim_rng", {
  m <- 4
  p <- 2
  q <- 3
  veclen <- 2 * p + m + q
  time <- 5
  iter <- 6
  x <- array(rep(seq_len(veclen), each = time * iter), c(iter, time, veclen))
  out <- ssm_extract(x, m, p, q, type = "sim_rng")
  expect_equal(out$y[1, 1, ], 1:2)
  expect_equal(out$a[1, 1, ], 3:6)
  expect_equal(out$eps[1, 1, ], 7:8)
  expect_equal(out$eta[1, 1, ], 9:11)
})

test_that("ssm_extract param option works", {
  m <- 4
  p <- 2
  # v, Finv, K, a, P
  veclen <- 1 + p + p * (p + 1) / 2 + m * p + m + m * (m + 1) / 2
  time <- 5
  iter <- 6
  x <- array(rep(seq_len(veclen), each = time * iter), c(iter, time, veclen))
  out <- ssm_extract(x, m, p, type = "filter", params = c("v", "K"))
  expect_length(out, 2)
  expect_named(out, c("v", "K"))
  expect_equal(out$v[1, 1, ], c(2, 3))
  expect_equal(out$K[1, 1, , ], matrix(6 + seq_len(m * p), m, p))
})