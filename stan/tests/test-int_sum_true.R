context("Stan function int_sum_true")

test_int_sum_true <- function(x, expected) {
  .data <- list(x = as.integer(x), n = length(x))
  modfit <- test_stan_function("int_sum_true", data = .data)
  output <- as.integer(rstan::extract(modfit, "output")[[1]])
  expect_equal(output, expected)
}

int_sum_true_params <-
  list(list(x = c(0L, 0L), expected = 0L),
       list(x = c(0L, 1L), expected = 1L),
       list(x = c(1L, 1L), expected = 2L),
       list(x = c(5L, 0L), expected = 1L),
       list(x = c(5L, 0L, -1L, 1L), expected = 2L)
  )

for (.x in int_sum_true_params) {
  test_that(sprintf("int_sum_true works with x = %s\n",
                    capture.output(dput(.x[["x"]]))), {
    invoke(test_int_sum_true, .x)
  })
}
