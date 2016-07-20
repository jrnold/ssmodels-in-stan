test_that("Stan file ssm.stan compiles and runs", {
  expect_equal(test_stan_function("ssm", output = FALSE), NULL)
})
