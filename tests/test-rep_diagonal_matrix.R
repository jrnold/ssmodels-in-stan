#function: rep_diagonal_matrix
context("rep_diagonal_matrix")

run_rep_diagonal_matrix <- function(x, m, n, k) {
  .data <- list(x = x, m = m, n = n, k = k)
  modfit <- test_stan_function("rep_diagonal_matrix", data = .data)
  rstan::extract(modfit, "output")[[1]][1, , ]
}

test_that("rep_diagonal_matrix works with x = 1, m = 2, n = 2, k = 0", {
  output <- run_rep_diagonal_matrix(1, 2, 2, 0)
  expected <- diag(1, 2, 2)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works x = 5, m = 2, n = 2, k = 0", {
  output <- run_rep_diagonal_matrix(5, 2, 2, 0)
  expected <- diag(5, 2, 2)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with  x = 5, m = 3, n = 2, k = 0", {
  output <- run_rep_diagonal_matrix(5, 3, 2, 0)
  expected <- diag(5, 3, 2)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with x = 5, m = 2, n = 3, k = 0", {
  output <- run_rep_diagonal_matrix(5, 2, 3, 0)
  expected <- diag(5, 2, 3)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with x = 2, m = 2, n = 3, k = 1", {
  output <- run_rep_diagonal_matrix(2, 2, 3, 1)
  expected <- matrix(c(0, 2, 0,
                       0, 0, 2), byrow = TRUE, 2, 3)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with x = 2, m = 2, n = 3, k = 5", {
  output <- run_rep_diagonal_matrix(2, 2, 3, 5)
  expected <- matrix(0, 2, 3)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with x = 2, m = 2, n = 3, k = -1", {
  output <- run_rep_diagonal_matrix(2, 2, 3, -1)
  expected <- matrix(c(0, 0, 0,
                       2, 0, 0), byrow = TRUE, 2, 3)
  expect_equal(output, expected)
})

test_that("rep_diagonal_matrix works with x = 2, m = 2, n = 3, k = -5", {
  output <- run_rep_diagonal_matrix(2, 2, 3, -5)
  expected <- matrix(0, 2, 3)
  expect_equal(output, expected)
})