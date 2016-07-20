test_that("Stan function select_indexes works", {
  f <- function(x) {
    x <- as.integer(x)
    m <- length(x)
    n <- sum(x > 0L)
    modfit <- test_stan_function("select_indexes",
                                 data = list(n = n, m = m, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(0, 1, 0, 3, -1)), as.integer(c(2, 4)))
  expect_equal(f(c(1, 1, 1)), as.integer(c(1, 2, 3)))

})
