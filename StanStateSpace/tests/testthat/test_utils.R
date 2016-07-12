context("test functions in utils.R")

test_that("vector to symmat produces symmetric matrices", {
  # 1 x 1
  expect_equal(structure(c(1L), .Dim = c(1L, 1L)),
               vector_to_symmat(1L))
  # 3 x 3
  expect_equal(vector_to_symmat(1:6),
               matrix(c(1, 2, 3,
                      2, 4, 5,
                      3, 5, 6), 3, 3))

})

test_that("vector_to_symmat produces matrices of the correct type", {
  expect_true(is.numeric(vector_to_symmat(c(1, 2, 3, 4))))
  expect_true(is.integer(vector_to_symmat(1:4)))
  expect_true(is.character(vector_to_symmat(c("a", "b", "c", "d"))))
  expect_true(is.logical(vector_to_symmat(c(TRUE, TRUE, FALSE, FALSE))))
})

test_that("as.Date.ts produces date objects", {
  x_months <- ts(start = 1900, end = 1901, frequency = 12)
  expect_equal(as.Date(x_months),
               as.Date(c("1900-01-01", "1900-02-01", "1900-03-01", "1900-04-01", "1900-05-01",
                        "1900-06-01", "1900-07-01", "1900-08-01", "1900-09-01", "1900-10-01",
                        "1900-11-01", "1900-12-01", "1901-01-01"),
                       format = "%Y-%m-%d"))
  x_qtrs <- ts(start = 1900, end = 1901, frequency = 4)
  expect_equal(as.Date(x_qtrs),
               as.Date(c("1900-01-01", "1900-04-01", "1900-07-01", "1900-10-01", "1901-01-01"),
               format = "%Y-%m-%d"))
  x_yrs <- ts(start = 1900, end = 1904, frequency = 1)
  expect_equal(as.Date(x_yrs),
               as.Date(c("1900-01-01", "1901-01-01", "1902-01-01", "1903-01-01", "1904-01-01"),
                         format = "%Y-%m-%d"))
  # 1900 not a leap year
  x_days <- as.Date(ts(start = 1900, end = 1901, frequency = 365L))[1:7]
  date_days <- as.Date(c("1900-01-01", "1900-01-01", "1900-01-02", "1900-01-04", "1900-01-05",
                         "1900-01-06", "1900-01-06"),
                       format = "%Y-%m-%d")
  expect_equivalent(x_days, date_days)
})

test_that("unit_lower_tri works", {
  expect_equal(unit_lower_tri(2),
               matrix(c(1, 1, 0, 1), 2, 2))
  expect_equal(unit_lower_tri(2, diag = FALSE),
               matrix(c(0, 1, 0, 0), 2, 2))
  expect_equal(unit_lower_tri(2, 3),
               matrix(c(1, 1, 0, 1, 0, 0), 2, 3))
})

test_that("unit_upper_tri works", {
  expect_equal(unit_upper_tri(2),
               matrix(c(1, 0, 1, 1), 2, 2))
  expect_equal(unit_upper_tri(2, diag = FALSE),
               matrix(c(0, 0, 1, 0), 2, 2))
  expect_equal(unit_upper_tri(2, 3),
               matrix(c(1, 0, 1, 1, 1, 1), 2, 3))
})

test_that("selection matrix works", {
  X <- matrix(1:9, 3, 3)
  expect_equal(selection_matrix(3, c(1, 3)) %*% X,
               X[c(1, 3), , drop = FALSE])
  expect_equal(selection_matrix(3, 1:3) %*% X, X)
  expect_equal(selection_matrix(3, 1) %*% X, X[1, , drop = FALSE])
  expect_equal(X %*% selection_matrix(3, c(1, 3), rows = FALSE),
               X[ , c(1, 3), drop = FALSE])
  expect_equal(X %*% selection_matrix(3, 1, rows = FALSE), X[ , 1, drop = FALSE])
  expect_equal(X %*% selection_matrix(3, 1:3, rows = FALSE), X)
})

test_that("eye works", {
  expect_equal(eye(2, 2), matrix(c(1, 0, 0, 1), 2, 2))
  expect_equal(eye(2, 2, -1), matrix(c(0, 1, 0, 0), 2, 2))
  expect_equal(eye(2, 2, 1), matrix(c(0, 0, 1, 0), 2, 2))
  expect_equal(eye(2, 2, 2), matrix(c(0, 0, 0, 0), 2, 2))
  expect_equal(eye(2, 2, -2), matrix(c(0, 0, 0, 0), 2, 2))
  expect_equal(eye(3, 2), matrix(c(1, 0, 0, 0, 1, 0), 3, 2))
  expect_equal(eye(3, 2, -1), matrix(c(0, 1, 0, 0, 0, 1), 3, 2))
  expect_equal(eye(3, 2, -2), matrix(c(0, 0, 1, 0, 0, 0), 3, 2))
  expect_equal(eye(3, 2, -3), matrix(0, 3, 2))
  expect_equal(eye(3, 2, 1), matrix(c(0, 0, 0, 1, 0, 0), 3, 2))
  expect_equal(eye(3, 2, 2), matrix(0, 3, 2))
})
