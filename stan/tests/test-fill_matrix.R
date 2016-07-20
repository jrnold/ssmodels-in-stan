test_that("Stan function fill_matrix works", {
  f <- function(x, m, n, i, j, a) {
    p <- nrow(x)
    q <- ncol(x)
    modfit <- test_stan_function("fill_matrix",
                                 data = list(m = m, n = n, p = p, q = q,
                                             x = x, i = array(as.integer(i)), j = array(as.integer(j)),
                                             a = a))
    ret <- rstan::extract(modfit)[["output"]]
    array(ret, dim(ret)[-1L])
  }
  m <- 4
  n <- 5
  x <- matrix(1:6, 2, 3)
  a <- 0.0
  i <- c(1, 3)
  j <- c(1, 2, 5)
  output <- f(x, m, n, i, j, a)
  expected <- structure(c(1, 0, 2, 0, 3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5,
                          0, 6, 0), .Dim = 4:5)
  expect_equal(output, expected)

})
