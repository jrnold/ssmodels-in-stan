context("Stan function matrix_pow")

test_matrix_pow <- function(m, A, n, expected) {
  .data <- list(m = m, A = A, n = n)
  modfit <- test_stan_function("matrix_pow", data = .data)
  output <- rstan::extract(modfit, "output")[[1]]
  dim(output) <- dim(output)[2:3]
  expect_equal(output, expected, TOL)
}

matrix_pow_params <- list()

matrix_pow_params[[1]] <- within(list(), {
  m <- 2
  A <- rand_mat(m, m)
  n <- 0
  expected <- diag(1, m, m)
  info <- "m = 2, n = 0"
})

matrix_pow_params[[2]] <- within(list(), {
  m <- 2
  A <- rand_mat(m, m)
  n <- 1
  expected <- A
  info <- "m = 2, n = 1"
})

matrix_pow_params[[3]] <- within(list(), {
  m <- 2
  A <- rand_mat(m, m)
  n <- 2
  expected <- A %*% A
  info <- "m = 2, n = 2"
})

matrix_pow_params[[4]] <- within(list(), {
  m <- 2
  A <- rand_mat(m, m)
  n <- 3
  expected <- A %*% A %*% A
  info <- "m = 2, n = 3"
})

matrix_pow_params[[4]] <- within(list(), {
  m <- 2
  A <- matrix(0, m, m)
  n <- 3
  expected <- A
  info <- "a zero matrix"
})

matrix_pow_params[[5]] <- within(list(), {
  m <- 5
  A <- rand_mat(3)
  n <- 3
  expected <- A %*% A %*% A
  info <- "m = 5, n = 3"
})

matrix_pow_params[[5]] <- within(list(), {
  m <- 5
  A <- diag(2, m, m)
  n <- 3
  expected <- diag(2 ^ n, m, m)
  info <- "diag matrix"
})

for (.x in matrix_pow_params) {
  test_that(paste0("matrix_pow works with ", .x$info), {
    test_matrix_pow(.x$m, .x$A, .x$n, .x$expected)
  })
}
