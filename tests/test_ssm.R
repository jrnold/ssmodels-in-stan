context("Stan models")

# Run a stan test
# @param file name of function. The file is "stan/test_{file}.stan"
# @param init Inputs to the function
# @param data Other data needed to initialize the model
# @return An array with the output of the function.
test_stan_function <- function(FUN, data, init = NULL) {
  filename <- test_path(file.path("stan",  paste0("test_", FUN, ".stan")))
  init <- if (!is.null(init)) list(init) else "random"
  ret <- ssm_stan_model(filename,
                        init = init, data = data, iter = 1, chains = 1,
                        algorithm = "Fixed_param")
  rstan::extract(ret, "output")[["output"]]
}

if (interactive()) {

  test_that("Stan function to_symmetric_matrix works", {
    f <- function(x) {
      ret <- test_stan_function("to_symmetric_matrix",
                         data = list(m = nrow(x), n = ncol(x), input = x))
      d <- dim(ret)
      array(ret, d[-1L])
    }
    expect_equal(f(matrix(1)), matrix(1))
    expect_equal(f(matrix(c(1, 2, 3, 4), 2, 2)),
                 matrix(c(1, 2.5, 2.5, 4), 2, 2))
  })

  test_that("Stan function to_matrix_colwise works", {
    f <- function(x, m, n) {
      ret <- test_stan_function("to_matrix_colwise",
                         data = list(m = m, n = n, input = x))
      array(ret, dim(ret)[-1L])
    }
    expect_equal(f(1:6, 2, 3), matrix(1:6, 2, 3))
  })

  test_that("Stan function symmat_size works", {
    f <- function(n) {
      ret <- test_stan_function("symmat_size", data = list(n = n), init = NULL)
      as.numeric(ret)
    }
    expect_equal(f(3), 6)
    expect_equal(f(2), 3)
    expect_equal(f(1), 1)
  })

  test_that("Stan function find_symmat_dim works", {
    f <- function(n) {
      ret <- test_stan_function("find_symmat_dim", data = list(n = n), init = NULL)
      as.numeric(ret)
    }
    expect_equal(f(1), 1)
    expect_equal(f(3), 2)
    expect_equal(f(6), 3)
  })

  test_that("Stan function vector_to_symmat works", {
    f <- function(x) {
      ret <- test_stan_function("vector_to_symmat",
                                data = list(x = x, m = length(x),
                                            n = floor(sqrt(2 * length(x)))),
                                init = NULL)
     ret
    }
    expect_equal(f(1:6), matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3))
    expect_equal(f(array(1)), matrix(1, 1, 1))
  })

  test_that("Stan function symmat_to_vector works", {
    f <- function(x) {
      ret <- test_stan_function("symmat_to_vector",
                                data = list(x = x,
                                            m = nrow(x),
                                            n = ncol(x)),
                                init = NULL)
      as.numeric(ret)
    }
    expect_equal(f(matrix(1, 1, 1)), 1)
    expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3)), as.numeric(1:6))
    # case with unequal rows/cols
    expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5), 3, 2)), c(1, 2, 4))
  })


  test_that("Stan function kronecker_prod works", {
    f <- function(A, B) {
      ret <- test_stan_function("kronecker_prod",
                                data = list(A = A, B = B,
                                            m = nrow(A), n = ncol(A),
                                            p = nrow(B), q = ncol(B)),
                                init = NULL)
      array(ret, dim(ret)[-1L])
    }
    A <- matrix(1:2, 1, 2)
    B <- matrix(3 + 1:12, 3, 4)
    expect_equal(f(A, B), kronecker(A, B))
  })

}
