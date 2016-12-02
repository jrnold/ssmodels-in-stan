#function: ssm_update_loglik_miss
context("ssm_update_loglik_miss")
# @param p number of observations
# @param missing logical vector with missing values indicated by true
test_ssm_update_loglik_miss <- function(p, missing) {
  v <- rnorm(p)
  Finv <- solve(rand_spd_mat(p))
  modfit <- test_stan_function("ssm_update_loglik_miss",
                               data = list(m = m, p = p, v = v, Finv = Finv))
  ret <- rstan::extract(modfit, "output")[[1]]
  expected <- -0.5 * p * log(2 * base::pi) - 0.5 * (log(det(solve(Finv))) + t(v) %*% Finv %*% v)
  output <- f(m, p, array(v), Finv)
  as.numeric(ret)
}

ssm_pdate_loglik_miss_param <-
  list(
    list(p = 1, missing = TRUE),
    list(p = 1, missing = FALSE),
    list(p = 2, missing = c(TRUE, FALSE)),
    list(p = 2, missing = c(TRUE, TRUE)),
    list(p = 2, missing = c(FALSE, FALSE)),
    list(p = 3, missing = c(FALSE, TRUE, FALSE))
    )

test_that("Stan function ssm_update_loglik_miss works", {

})
