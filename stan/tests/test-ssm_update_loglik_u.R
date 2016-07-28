#function: ssm_update_loglik_u
context("ssm_update_loglik")

test_ssm_update_loglik_u <- function(v, Finv) {
  modfit <- test_stan_function("ssm_update_loglik_u",
                               data = list(v = v, Finv = Finv))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- -0.5 * (log(2 * base::pi) - log(Finv) + Finv * v ^ 2)
  expect_equal(output, as.numeric(expected), tol = TOL)
}

ssm_update_loglik_u_params <-
  list(
      list(v = rnorm(1), Finv = exp(rnorm(1)))
  )

for (.x in ssm_update_loglik_u_params) {
  msg <- sprintf("It works with v = %f, Finv = %f", .x$v, .x$Finv)
  test_that(msg, {
    test_ssm_update_loglik_u(.x$v, .x$Finv)
  })
}
