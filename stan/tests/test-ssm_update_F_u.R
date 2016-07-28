#function: ssm_update_F_u
context("ssm_update_F_u")

test_ssm_update_Finv_u <- function(m, Z, P, H) {
  modfit <-
    test_stan_function("ssm_update_F_u",
                       data = list(m = m, Z = Z, P = P, H = H))
  output <- as.numeric(rstan::extract(modfit, "output")[[1]])
  expected <- as.numeric(t(Z) %*% P %*% Z + H)
  expect_equal(output, expected, TOL)
}

rand_ssm_update_Finv_u_params <- function(m) {
  list(m = m, Z = array(rnorm(m)),
       P = rand_spd_mat(m), H = exp(rnorm(1)))
}

test_ssm_update_Finv_u_params <- map(c(1, 2, 3), rand_ssm_update_Finv_u_params)

for (.x in test_ssm_update_Finv_u_params) {
  test_ssm_update_Finv_u(.x$m, .x$Z, .x$P, .x$H)
}
