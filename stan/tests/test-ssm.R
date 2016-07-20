context("Stan models")

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

test_that("Stan function symmat_to_vector works", {
  f <- function(x) {
    modfit <- test_stan_function("symmat_to_vector",
                              data = list(x = x,
                                          m = nrow(x),
                                          n = ncol(x)),
                              init = NULL)
    as.numeric(rstan::extract(modfit, "output")[[1]])
  }
  expect_equal(f(matrix(1, 1, 1)), 1)
  expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5, 3, 5, 6), 3, 3)), as.numeric(1:6))
  # case with unequal rows/cols
  expect_equal(f(matrix(c(1, 2, 3, 2, 4, 5), 3, 2)), c(1, 2, 4))
})


test_that("Stan function kronecker_prod works", {
  f <- function(A, B) {
    modfit <- test_stan_function("kronecker_prod",
                              data = list(A = A, B = B,
                                          m = nrow(A), n = ncol(A),
                                          p = nrow(B), q = ncol(B)),
                              init = NULL)
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  A <- matrix(1:2, 1, 2)
  B <- matrix(3 + 1:12, 3, 4)
  expect_equal(f(A, B), kronecker(A, B))
})


#' Update a should calculate T * a + K * v + c
#' The output should be m x 1
test_that("Stan function ssm_filter_update_a works", {
  f <- function(m, p, a, c, T, v, K) {
    modfit <- test_stan_function("ssm_filter_update_a",
                                 data = list(m = m, p = p, a = a,
                                             c = c, T = T, v = v, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  a <- rnorm(m)
  c <- rnorm(m)
  T <- matrix(rnorm(m * m), m, m)
  v <- rnorm(p)
  K <- matrix(runif(m * p), m, p)
  expected <- T %*% a + K %*% v + c
  anew <- f(m, p, array(a), array(c), T, array(v), K)
  expect_length(anew, m)
  expect_equal(matrix(anew), expected, tolerance = 10e-5)
})

#' Update a should calculate (DK, p. 107)
#' $$
#' P_{t + 1} = T P_t (T - K Z)' + R Q R
#' $$
test_that("Stan function ssm_filter_update_P works", {
  f <- function(m, p, P, Z, T, RQR, K) {
    modfit <- test_stan_function("ssm_filter_update_P",
                                 data = list(m = m, p = p, P = P,
                                             Z = Z, T = T, RQR = RQR, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  P <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  T <- matrix(rnorm(m * m), m, m)
  # Need Kalman gain to be generated in this way; otherwise (T - K Z) not symm
  K <- T %*% P %*% t(Z) %*% solve(rand_pdmat(p))
  RQR <- rand_pdmat(m)
  expected <- T %*% P %*% t(T - K %*% Z) + RQR
  output <- f(m, p, P, Z, T, RQR, K)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})

#' Update a should calculate (DK, Sec 4.3.2)
#' $$
#' F = Z * P * Z' + H
#' $$
test_that("Stan function ssm_filter_update_F works", {
  f <- function(m, p, Z, P, H) {
    modfit <- test_stan_function("ssm_filter_update_F",
                                 data = list(m = m, p = p, Z = Z, P = P, H = H))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_pdmat(m)
  H <- rand_pdmat(p)
  expected <- Z %*% P %*% t(Z) + H
  output <- f(m, p, Z, P, H)
  expect_length(output, p * p)
  expect_equal(output, expected, tolerance = 10e-5)
})


#' #' Update a should calculate (DK, Sec 4.3.2)
#' $$
#' F^{-1} = (Z * P * Z' + H)^{-1}
#' $$
test_that("Stan function ssm_filter_update_Finv works", {
  f <- function(m, p, Z, P, H) {
    modfit <- test_stan_function("ssm_filter_update_Finv",
                                 data = list(m = m, p = p, Z = Z, P = P, H = H))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_pdmat(m)
  H <- rand_pdmat(p)
  expected <- solve(Z %*% P %*% t(Z) + H)
  output <- f(m, p, Z, P, H)
  expect_length(output, p * p)
  expect_equal(output, expected, tolerance = 10e-5)
})

test_that("Stan function ssm_filter_update_Finv_miss works", {
  f <- function(m, p, p_t, Z, P, H, y_idx) {
    modfit <- test_stan_function("ssm_filter_update_Finv_miss",
                                 data = list(m = m, p = p, p_t = p_t,
                                   Z = Z, P = P,
                                   H = H, y_idx = y_idx))
    if (p_t == 0) {
      expected <- matrix(0, p, p)
    } else if (p_t == p) {
      expected <- solve(Z %*% P %*% t(Z) + H)
    } else {
      idx <- y_idx[1:p_t]
      Z_star <- Z[idx, , drop = FALSE]
      H_star <- H[idx, idx, drop = FALSE]
      expected <- matrix(0, p, p)
      Finv_star <- solve(Z_star %*% P %*% t(Z_star) + H_star)
      for (i in seq_along(idx)) {
        for (j in seq_along(idx)) {
          expected[idx[i], idx[j]] <- Finv_star[i, j]
        }
      }
    }
    output <- rstan::extract(modfit, "output")[[1]][1, , ]
    expect_length(output, p * p)
    expect_equal(output, expected, tolerance = 10e-5)
  }

  # No missing values
  m <- 3L
  p <- 4L
  p_t <- p
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_pdmat(m)
  H <- rand_pdmat(p)
  f(m, p, p_t, Z, P, H, y_idx = seq_len(p_t))

  # Some missing values
  p_t <- 2
  f(m, p, p_t, Z, P, H, y_idx = c(1, 3, 0, 0))

  # All missing values
  p_t <- 0
  f(m, p, p_t, Z, P, H, y_idx = rep(0, p))

})

test_that("Stan function ssm_filter_update_v works", {
  f <- function(m, p, y, a, d, Z) {
    modfit <- test_stan_function("ssm_filter_update_v",
                                 data = list(m = m, p = p,
                                   a = array(a), d = array(d), y = array(y),
                                   Z = Z))

    expected <- as.numeric(y - d - Z %*% matrix(a))
    output <- as.numeric(rstan::extract(modfit, "output")[[1]][1, ])
    expect_length(output, p)
    expect_equal(output, expected, tolerance = 10e-5)
  }
  p <- 4L
  m <- 3L
  y <- rnorm(p)
  a <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  d <- rnorm(p)
  f(m, p, y, a, d, Z)
})

test_that("Stan function ssm_filter_update_v_miss works", {
  f <- function(m, p, y, a, d, Z, p_t, y_idx) {
    modfit <- test_stan_function("ssm_filter_update_v_miss",
                                 data = list(m = m, p = p, p_t = p_t,
                                   a = array(a),  d = array(d), y = array(y),
                                   Z = Z, y_idx = y_idx))
    expected <- rep(0, p)
    for (i in seq_len(p)) {
      if (i %in% y_idx) {
        expected[i] <- as.numeric(y[i] - d[i] - Z[i, ] %*% matrix(a))
      }
    }
    output <- as.numeric(rstan::extract(modfit, "output")[[1]][1, ])
    print(output)
    print(expected)
    expect_length(output, p)
    expect_equal(output, expected, tolerance = 10e-5)
  }

  # No missing values
  m <- 3L
  p <- 4L
  p_t <- p
  y <- rnorm(p)
  a <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  d <- rnorm(p)
  f(m, p, y, a, d, Z, p_t, y_idx = seq_len(p))

  # Some missing values
  p_t <- 2
  f(m, p, y, a, d, Z, p_t, y_idx = c(1, 3, 0, 0))

  # All missing values
  p_t <- 0
  f(m, p, y, a, d, Z, p_t, y_idx = rep(0, p))

})

#' Update of K should calculate (DK, Sec 4.3.2)
#' $$
#' K = T * P * Z' * F^{-1}
#' $$
test_that("Stan function ssm_filter_update_K works", {
  f <- function(m, p, P, Z, T, Finv) {
    modfit <- test_stan_function("ssm_filter_update_K",
                                 data = list(m = m, p = p, T = T, Z = Z,
                                             P = P, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  T <- matrix(rnorm(m * m), m, m)
  Z <- matrix(rnorm(m * p), p, m)
  P <- rand_pdmat(m)
  Finv <- solve(rand_pdmat(p))
  expected <- T %*% P %*% t(Z) %*% Finv
  output <- f(m, p, P, Z, T, Finv)
  expect_length(output, m * p)
  expect_equal(output, expected, tolerance = 10e-5)
})

#' Update of L should calculate (DK, Sec 4.3.5)
#' $$
#' L = T - K * Z
#' $$
test_that("Stan function ssm_filter_update_L works", {
  f <- function(m, p, Z, T, K) {
    modfit <- test_stan_function("ssm_filter_update_L",
                                 data = list(m = m, p = p, Z = Z, T = T, K = K))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  m <- 3L
  Z <- matrix(rnorm(m * p), p, m)
  T <- matrix(rnorm(m * m), m, m)
  P <- rand_pdmat(m)
  # Need Kalman gain to be generated in this way; otherwise (T - K Z) not symm
  K <- T %*% P %*% t(Z) %*% solve(rand_pdmat(p))
  expected <- T - K %*% Z
  output <- f(m, p, Z, T, K)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})

#' Update of log-likelihood is
#' $$
#' L_t = - 1/p log(2 * pi) - 1/2 (log(|F|) + v' * F^{-1} * v)
#' $$
test_that("Stan function ssm_filter_update_ll works", {
  f <- function(m, p, v, Finv) {
    modfit <- test_stan_function("ssm_filter_update_ll",
                                 data = list(m = m, p = p, v = v, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    as.numeric(ret)
  }
  m <- 3L
  p <- 4L
  m <- 3L
  v <- rnorm(p)
  Finv <- solve(rand_pdmat(p))
  expected <- -0.5 * p * log(2 * base::pi) - 0.5 * (log(det(solve(Finv))) + t(v) %*% Finv %*% v)
  output <- f(m, p, array(v), Finv)
  expect_length(output, 1)
  expect_equal(output, as.numeric(expected), tol = 10e-5)
})

test_that("Stan function ssm_filter_idx and ssm_filter_get_* functions work", {


  f <- function(m, p, x) {
    modfit <- test_stan_function("ssm_filter_get",
                                 data = list(m = m, p = p, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  symmat_size <- function(n) n * (n + 1) / 2
  xsz <- 1 + p + symmat_size(p) + m * p + m + symmat_size(m)
  x <- seq_len(xsz)
  output <- f(m, p, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(output[["idx"]][1, , ],
               matrix(c(1, 1, 1,
                      p, 1 + 1, 1 + p,
                      symmat_size(p), 2 + p, 1 + p + symmat_size(p),
                      m * p, 2 + p + symmat_size(p), 1 + p + symmat_size(p) + m * p,
                      m, 2 + p + symmat_size(p) + m * p, 1 + p + symmat_size(p) + m * p + m,
                      symmat_size(m), 2 + p + symmat_size(p) + m * p + m, 1 + p + symmat_size(p) + m * p + m + symmat_size(m)),
                      byrow = TRUE, 6, 3))
  expect_equivalent(output[["ll"]], array(1))
  expect_equivalent(output[["v"]][1, ], 2:5)
  expect_equivalent(output[["Finv"]][1, , ],
                    matrix(c(6, 7, 8, 9,
                             7, 10, 11, 12,
                             8, 11, 13, 14,
                             9, 12, 14, 15),
                           4, 4))
  expect_equivalent(output[["K"]][1, ,],
                    matrix(16:(16 + m * p - 1), m, p))
  expect_equivalent(as.numeric(output[["a"]]),
                    27 + 1:m)
  expect_equivalent(output[["P"]][1, ,],
                    matrix(c(31, 32, 33,
                             32, 34, 35,
                             33, 35, 36),
                           3, 3))
})

#' Update a should calculate
#' $$
#' \E(\vec{alpha}_t | \vec{y}_{1:t}) = \vec{a}_{t|t} = \mat{T}_t * \vec{a}_t + \mat{K}_t \vec{v}_t .
#' $$
#' The output should be m x 1
test_that("Stan function ssm_filter_states_update_a works", {
  f <- function(m, p, a, P = P, Z, v, Finv) {
    modfit <- test_stan_function("ssm_filter_states_update_a",
                                 data = list(m = m, p = p, a = a, P = P,
                                             Z = Z, v = v, Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  a <- rnorm(m)
  P <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  v <- rnorm(p)
  Finv <- solve(rand_pdmat(p))
  expected <- a + P %*% t(Z) %*% Finv %*% v
  output <- f(m, p, a, P, Z, v, Finv)
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})

#' Update should calculate (DK, Sec 4.3.2)
#' $$
#' \Var(\vec{alpha}_t | \vec{y}_{1:t}) = \mat{P}_{t|t} = \mat{P}_t - \mat{P}_t \mat{Z}_t' \mat{F}_t^{-1} \mat{Z}_t \mat{P}_t .
#' $$
test_that("Stan function ssm_filter_states_update_P works", {
  f <- function(m, p, P, Z, Finv) {
    modfit <- test_stan_function("ssm_filter_states_update_P",
                                 data = list(m = m, p = p, P = P, Z = Z,
                                             Finv = Finv))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  P <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  Finv <- solve(rand_pdmat(p))
  expected <- P - P %*% t(Z) %*% Finv %*% Z %*% P
  output <- f(m, p, P, Z, Finv)
  expect_length(output, m * m)
  expect_equal(output, expected, tolerance = 10e-5)
})

test_that("Stan functions ssm_filter_states_get", {
  f <- function(m, x) {
    modfit <- test_stan_function("ssm_filter_states_get",
                                 data = list(m = m, p = p, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  symmat_size <- function(n) n * (n + 1) / 2
  xsz <- m + symmat_size(3)
  x <- seq_len(xsz)
  output <- f(m, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["a"]]), 1:3)
  expect_equal(output[["P"]][1, , ],
               matrix(c(4, 5, 6,
                        5, 7, 8,
                        6, 8, 9), 3, 3))
})

test_that("Stan function matrix_diff works", {
  f <- function(m, n, A, B) {
    modfit <- test_stan_function("matrix_diff",
                                 data = list(m = m, n = n,
                                             A = A, B = B))
    rstan::extract(modfit)[["output"]]
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  n <- 4L
  A <- matrix(rnorm(m * n), m, n)
  B <- matrix(rnorm(m * n), m, n)
  output <- f(m, n, A, B)
  expected <- max(abs(A - B)) / max(abs(A))
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})

#' Expected value is in DK Sec 4.4.4
test_that("Stan function ssm_smoth_update_r works", {
  f <- function(m, p, r, Z, v, Finv, L) {
    modfit <- test_stan_function("ssm_smooth_update_r",
                                 data = list(m = m, p = p, r = r, Z = Z,
                                             v = v, Finv = Finv, L = L))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  r <- rnorm(m)
  Z <- matrix(rnorm(m * p), p, m)
  v <- rnorm(p)
  Finv <- solve(rand_pdmat(p))
  L <- matrix(rnorm(m * m), m, m)
  expected <- t(Z) %*% Finv %*% v + t(L) %*% r
  output <- f(m, p, array(r), Z, array(v), Finv, L)
  expect_equal(as.numeric(output), as.numeric(expected), tol = 10e-5)
})

#' Expected value is in DK Sec 4.4.4
test_that("Stan function ssm_smoth_update_N works", {
  f <- function(m, p, N, Z, Finv, L) {
    modfit <- test_stan_function("ssm_smooth_update_N",
                                 data = list(m = m, p = p,
                                             N = N,
                                             Z = Z,
                                             Finv = Finv, L = L))
    ret <- rstan::extract(modfit, "output")[[1]]
    array(ret, dim(ret)[-1L])
  }
  m <- 3L
  p <- 4L
  N <- rand_pdmat(m)
  Z <- matrix(rnorm(m * p), p, m)
  Finv <- solve(rand_pdmat(p))
  L <- matrix(rnorm(m * m), m, m)
  expected <- t(Z) %*% Finv %*% Z + t(L) %*% N %*% L
  output <- f(m, p, N, Z, Finv, L)
  expect_equal(output, expected, tol = 10e-5)
})

test_that("Stan functions ssm_filter_states_get", {
  f <- function(m, x) {
    modfit <- test_stan_function("ssm_smooth_state_get",
                                 data = list(m = m, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  xsz <- 9
  x <- seq_len(xsz)
  output <- f(m, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["alpha"]]), 1:3)
  expect_equal(output[["V"]][1, , ],
               matrix(c(4, 5, 6,
                        5, 7, 8,
                        6, 8, 9), 3, 3))
})

test_that("Stan functions ssm_filter_eps_get works", {
  f <- function(p, x) {
    modfit <- test_stan_function("ssm_smooth_eps_get",
                                 data = list(p = p, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 3L
  p <- 4L
  xsz <- 14
  x <- seq_len(xsz)
  output <- f(p, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["eps_mean"]]), 1:4)
  expect_equal(output[["eps_var"]][1, , ],
               matrix(c(5, 6,  7,  8,
                        6, 9,  10, 11,
                        7, 10, 12, 13,
                        8, 11, 13, 14), 4, 4, byrow = TRUE))
})

test_that("Stan functions ssm_filter_eta_get works", {
  f <- function(q, x) {
    modfit <- test_stan_function("ssm_smooth_eta_get",
                                 data = list(q = q, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  q <- 3L
  p <- 4L
  xsz <- 9
  x <- seq_len(xsz)
  output <- f(q, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(as.numeric(output[["eta_mean"]]), 1:3)
  expect_equal(output[["eta_var"]][1, , ],
               matrix(c(4, 5,  6,
                        5, 7,  8,
                        6, 8, 9), 3, 3, byrow = TRUE))
})

test_that("Stan functions ssm_sim_get works", {
  f <- function(m, p, q, x) {
    modfit <- test_stan_function("ssm_sim_get",
                                 data = list(m = m, p = p, q = q, x = x,
                                             xsz = length(x)))
    rstan::extract(modfit)
  }
  # m = 3, p = 4 are hardcoded in some later tests
  m <- 5L;
  p <- 4L;
  q <- 3L;
  xsz <- 2 * p + m + q
  x <- seq_len(xsz)
  output <- f(m, p, q, x)
  expect_equal(as.numeric(output[["sz"]]), length(x))
  expect_equal(output[["idx"]][1, , ],
               matrix(c(4, 1, 4,
                        5, 5, 9,
                        4, 10, 13,
                        3, 14, 16), 4, 3, byrow = TRUE))
  expect_equal(as.numeric(output[["y"]]), 1:4)
  expect_equal(as.numeric(output[["a"]]), 5:9)
  expect_equal(as.numeric(output[["eps"]]), 10:13)
  expect_equal(as.numeric(output[["eta"]]), 14:16)

})

#' Compare results of Stationary coefficients against R's implementation of the transformation

#' Transform vector of real numbers to stationary AR(p) coefficients
#'
#' Extracted from `stats::arima()`
ar_trans <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI). I just need to worry about one set of coefficients
  # so ignore the other parts
  .Call(stats:::C_ARIMA_transPars, as.numeric(par), as.integer(c(length(par), rep(0, 5))), TRUE)[[1]]
}

#' Transform PACF to ACF
#'
#' Extracted from `stats::arima()`
pacf_to_acf <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI)
  # Undo the tanh transformation in C_ARIMA_tranPars to get back to partial autocorrelation
  ar_trans(atanh(par))
}

#' Transform vector of stationary AR(p) coefficients to real numbers
#'
#' Extracted from `stats::arima()`
ar_invtrans <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI)
  # goes from AR to pacf, pacf (-1, 1) to (-infty, infty)
  .Call(stats:::C_ARIMA_Invtrans, as.numeric(par), as.integer(c(length(par), rep(0, 5))))
}

#' Transform ACF to PACF
acf_to_pacf <- function(par) {
  # Redo the tanh transformation to get to partial autocorrelations
  tanh(ar_invtrans(par))
}

#' Stationary ARMA covariance
#'
#' Extracted parts from the function `stats::arima`.
arma_init <- function(theta, phi, method = "Gardner1980", tol = 0) {
  if (method == "Gardner1980") {
    Q0 <- .Call(stats:::C_getQ0, as.numeric(phi), as.numeric(theta))
  } else {
    Q0 <- .Call(stats:::C_getQ0bis, as.numeric(phi), as.numeric(theta), tol = tol)
  }
  Q0
}

test_that("Stan function constrain_stationary works", {
  f <- function(x) {
    modfit <- test_stan_function("constrain_stationary", data = list(n = length(x), x = array(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- rnorm(i)
    output <- f(x)
    expected <- ar_trans(x)
    expect_equal(output, expected, tolerance = 10e-5)
  }
})

test_that("Stan function pacf_to_acf works", {
  f <- function(x) {
    modfit <- test_stan_function("pacf_to_acf", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- tanh(rnorm(i))
    expected <- pacf_to_acf(x)
    output <- f(x)
    expect_equal(output, expected, tolerance = TOL)
  }
})

test_that("Stan function unconstrain_stationary works", {
  f <- function(x) {
    modfit <- test_stan_function("unconstrain_stationary", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    # real numbers
    expected <- rnorm(i)
    # AR coefficients
    phi <- ar_trans(expected)
    output <- f(phi)
    expect_equal(f(phi), expected, tolerance = TOL)
  }
})

test_that("Stan function acf_to_pacf works", {
  f <- function(x) {
    modfit <- test_stan_function("acf_to_pacf", data = list(x = array(x), n = length(x)))
    as.numeric(rstan::extract(modfit)[["output"]])
  }
  for (i in 1:3) {
    x <- ar_trans(rnorm(i))
    expect_equal(f(x), acf_to_pacf(x), tolerance = TOL)
  }
})

test_that("Stan function ssm_stationary_cov works", {
  f <- function(T, RQR) {
    modfit <- test_stan_function("stationary_cov",
                                 data = list(T = T, RQR = RQR,
                                             m = nrow(T)))
    ret <- rstan::extract(modfit)[["output"]]
    array(ret, dim(ret)[-1L])
  }

  # ARMA(1, 1)
  arma11 <- function(phi, theta) {
    x <- matrix(0, 2, 2)
    x[1, 1] <- (1 + theta ^ 2 + 2 * theta * phi) / (1 - phi ^  2)
    x[2, 1] <- x[1, 2] <- theta
    x[2, 2] <- theta ^ 2
    x
  }
  phi <- 0.5
  theta <- 0.25
  arma11(0.5, 0.25)
  T <- matrix(c(phi, 0, 1, 0), 2, 2)
  R <- matrix(c(1, theta))
  RQR <- tcrossprod(R)
  expect_equal(f(T, RQR), arma11(phi, theta))

  # AR(1, 1)
  phi <- 0.5
  T <- matrix(phi)
  RQR <- matrix(1)
  f(T, RQR)
  expect_equal(as.numeric(f(T, RQR)), 1 / (1 - phi ^ 2), tolerance = TOL)

})

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

test_that("Stan function sum_true works", {
  f <- function(x) {
    n <- length(x)
    modfit <- test_stan_function("int_sum_true",
                                 data = list(n = n, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(1, -1, 3, 0, 0)), 2L)
  expect_equal(f(c(1, 1, 1)), 3L)
  expect_equal(f(c(0, 0, 0)), 0L)

})

test_that("Stan function sum_false works", {
  f <- function(x) {
    n <- length(x)
    modfit <- test_stan_function("int_sum_false",
                                 data = list(n = n, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(1, -1, 3, 0, 0)), 3L)
  expect_equal(f(c(1, 1, 1)), 0L)
  expect_equal(f(c(0, 0, 0)), 3L)

})

test_that("Stan function mask_indexes works", {
  f <- function(x) {
    x <- as.integer(x)
    m <- length(x)
    n <- sum(x <= 0L)
    modfit <- test_stan_function("mask_indexes",
                                 data = list(n = n, m = m, x = array(x)))
    ret <- rstan::extract(modfit)[["output"]]
    as.integer(ret)
  }
  expect_equal(f(c(1, -1, 3, 0, 0)), as.integer(c(2, 4, 5)))
  expect_equal(f(c(0, 0, 0)), as.integer(c(1, 2, 3)))

})

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


test_that("Stan function cholesky_decompose2 works", {
  f <- function(A) {
    n <- nrow(A)
    modfit <- test_stan_function("cholesky_decompose2",
                                 data = list(n = n, A = A))
    rstan::extract(modfit)[["output"]][1, , ]
  }
  A <- rand_pdmat(2)
  cholA <- t(chol(A))
  expect_equal(f(A), cholA, tolerance = 1e-5)
  A1 <- A[1, 1]
  cholA1 <- cholA[1, 1]
  expect_equal(f(matrix(c(A1, 0, 0, 0), 2, 2)),
                 matrix(c(cholA[1, 1], 0, 0, 0), 2, 2), tolerance = 1e-5)
  expect_equal(f(matrix(0, 2, 2)), matrix(0, 2, 2), tolerance = 1e-5)
})

test_that("Stan function normal2_rng works", {
  f <- function(n, mu, sigma) {
    modfit <- test_stan_function("normal2_rng",
                                 data = list(n = n, mu = mu, sigma = sigma))
    output <- as.numeric(rstan::extract(modfit)[["output"]])
    if (sigma == 0) {
      expect_true(all(output == mu))
    } else {
      expect_true(sd(output) > 0)
      z <- (mu - mean(output)) / (sigma / sqrt(n))
      expect_true(abs(z) < 4)
      expect_true(abs(sd(output) / sigma / sqrt(2 * (n - 1))) < 4)
    }
  }
  f(1000, 1, 2)
  f(10, 1, 0)
})

test_that("Stan function multi_normal2_rng works", {
  f <- function(n, mu, Sigma) {
    m <- length(mu)
    modfit <- test_stan_function("multi_normal2_rng",
                                 data = list(n = n, m = m,
                                             mu = mu, Sigma = Sigma))
    output <- rstan::extract(modfit)[["output"]][1, , ]
    print(dim(output))
    Sigma_zero <- (Sigma == 0)
    for (i in 1:m) {
      if (Sigma_zero[i, i]) {
        expect_true(all(output[, i] == mu[i]))
      } else {
        expect_true(sd(output[, i]) > 0)
        z <- (mu[i] - mean(output[ , i])) / (sqrt(Sigma[i, i] / sqrt(n)))
        print(z)
        expect_true(abs(z) < 4)
      }
    }
  }
  f(1000, c(1, 2), rand_pdmat(2))
  f(1000, c(1, 2), matrix(c(1, 0, 0, 0), 2, 2))
  f(1000, c(1, 2), matrix(0, 2, 2))
})

test_that("Stan function multi_normal_cholesky2_rng works", {
  f <- function(n, mu, L) {
    m <- length(mu)
    modfit <- test_stan_function("multi_normal_cholesky2_rng",
                                 data = list(n = n, m = m,
                                             mu = mu, L = L))
    output <- rstan::extract(modfit)[["output"]][1, , ]
    Sigma <- L %*% t(L)
    Sigma_zero <- (Sigma == 0)
    for (i in 1:m) {
      if (Sigma_zero[i, i]) {
        expect_true(all(output[, i] == mu[i]))
      } else {
        expect_true(sd(output[, i]) > 0)
        z <- (mu[i] - mean(output[ , i])) / (sqrt(Sigma[i, i] / sqrt(n)))
        print(z)
        expect_true(abs(z) < 4)
      }
    }
  }
  f(1000, c(1, 2), t(chol(rand_pdmat(2))))
  f(1000, c(1, 2), matrix(c(1, 0, 0, 0), 2, 2))
  f(1000, c(1, 2), matrix(0, 2, 2))
})
