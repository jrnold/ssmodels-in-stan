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
