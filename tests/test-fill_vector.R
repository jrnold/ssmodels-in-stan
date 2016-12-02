#function: fill_vector
context("fill_vector")
test_that("Stan function fill_vector works", {
  f <- function(x, n, i, a) {
    m <- length(x)
    modfit <- test_stan_function("fill_vector",
                                 data = list(m = m, n = n,
                                             x = x, i = array(as.integer(i)),
                                             a = a))
    ret <- rstan::extract(modfit)[["output"]]
    as.numeric(ret)
  }
  n <- 5
  x <- 1:3
  a <- 0.0
  i <- c(1, 3, 4)
  output <- f(x, n, i, a)
  expected <- c(1, 0, 2, 3, 0)
  expect_equal(output, expected)
})
