context("Stan loglikelihood tests")
set.seed(60322390)

test_that("Stan function ssm_loglik works", {
  f <- function(y, .d, .Z, .H, .c, .T, .R, .Q, a1, P1) {
    .data <- list(y = y,
                  n = nrow(y),
                  p = ncol(y),
                  m = ncol(.Z),
                  q = ncol(.Q),
                  d = .d, Z = .Z, H = .H, c = .c, T = .T, R = .R, Q = .Q,
                  a1 = a1, P1 = P1)
    modfit <- test_stan_function("ssm_loglik_constant", data = .data)
    output <- rstan::extract(modfit, "output")[[1]]
    ssm_extract(output, m = ncol(.Z), p = ncol(y), type = "filter")
  }

  sigma2_alpha <- 1469
  sigma2_eps <- 15099

  data("Nile", package = "datasets")
  output <- f(matrix(as.numeric(Nile)),
    .d = array(0, dim = c(1, 1)),
    .Z = array(1, dim = c(1, 1, 1)),
    .H = array(sigma2_alpha, dim = c(1, 1, 1)),
    .c = array(0, dim = c(1, 1)),
    .T = array(1, dim = c(1, 1, 1)),
    .R = array(1, dim = c(1, 1, 1)),
    .Q = array(sigma2_alpha, dim = c(1, 1, 1)),
    a1 = array(0),
    P1 = matrix(1e6))


})

