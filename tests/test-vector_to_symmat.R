#function: vector_to_symmat
context("vector_to_symmat")
test_that("Stan function vector_to_symmat works", {
  f <- function(x) {
    modfit <- test_stan_function("vector_to_symmat",
                              data = list(x = x, m = length(x),
                                          n = floor(sqrt(2 * length(x)))),
                              init = NULL)
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1])
  }
  expect_equal(f(1:6), matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3))
  expect_equal(f(array(1)), matrix(1, 1, 1))
})
