context("test ssm functions")

test_that("ssm_extract throws errors for invalid p, q, m arguments", {
  p <- 2
  m <- 4
  q <- 3
  iter <- 3
  time <- 2
  veclen <- 27
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
               regexp = "dim\\(x\\)\\[3\\] not equal to")
  expect_error(ssm_extract(x, m = m, p = p, q = q,
                           params = c("zzz")),
               regexp = "‘param’ value.*are invalid")

})

test_that("ssm_extract works for ...", {

})
