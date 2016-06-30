library("rstan")
library("testthat")
library("rprojroot")

expect_rstan_model_parses <- function(filename, ...) {
  include_dir <- find_rstudio_root_file("stan", "includes")
  exp <- capture_output(stanc_builder(file = filename, isystem = include_dir))
}


for (stan_file in dir(find_rstudio_root_file("tests", "models"),
                      pattern = "\\.stan$",
                      full.names = TRUE)) {
  print(stan_file)
  test_that(paste0(stan_file, " compiles"), {
    x <- run_stan_model(stan_file)
    expect_false(is(x$r, "try-error"))
  })
}


