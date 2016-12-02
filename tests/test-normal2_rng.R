#function: normal2_rng
context("normal2_rng")
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
