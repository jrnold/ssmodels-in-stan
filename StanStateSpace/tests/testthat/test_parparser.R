context("MCMC Paramter parsers")


test_that("expand_grid_dim works", {
  output <- StanStateSpace:::expand_grid_dim(c(2, 3, 4))
  expected <- data_frame(dim_1 = rep(1:2, times = 12),
                         dim_2 = rep(1:3, each = 2, times = 4),
                         dim_3 = rep(1:4, each = 6, times = 1))
  expect_equivalent(output, expected)
})

test_that("expand_grid_dim works with colwise = FALSE", {
  output <- StanStateSpace:::expand_grid_dim(c(2, 3, 4), colmajor = FALSE)
  expected <- data_frame(dim_1 = rep(1:2, each = 12, times = 1),
                         dim_2 = rep(1:3, each = 4, times = 2),
                         dim_3 = rep(1:4, times = 6))
  expect_equivalent(output, expected)
})

test_that("create_stan_parnames works as expected", {
  expect_equal(create_stan_parnames("alpha", NA), "alpha")
  expect_equal(create_stan_parnames("alpha", NULL), "alpha")
  expect_equal(create_stan_parnames("alpha", 0L), "alpha")
  expect_equal(create_stan_parnames("alpha", 1L), "alpha[1]")
  expect_equal(create_stan_parnames("alpha", 2L), c("alpha[1]", "alpha[2]"))
  expect_equal(create_stan_parnames("alpha", c(2L, 2L)),
               c("alpha[1,1]", "alpha[2,1]", "alpha[1,2]", "alpha[2,2]"))
})

test_that("create_parnames works as expected", {
  expect_equal(create_parnames("alpha", 0L, "<", ";", ">"), "alpha")
  expect_equal(create_parnames("alpha", 1L, "<", ";", ">"), "alpha<1>")
  expect_equal(create_parnames("alpha", 2L, "<", ";", ">"),
               c("alpha<1>", "alpha<2>"))
  expect_equal(create_parnames("alpha", c(2L, 2L), "<", ";", ">"),
               c("alpha<1;1>", "alpha<2;1>", "alpha<1;2>", "alpha<2;2>"))
})
