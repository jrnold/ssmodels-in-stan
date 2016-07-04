context("Stan models")

if (interactive()) {
  test_that("ssm.stan parses and compiles", {
    isystem <- system.file("stan/include", package = "StanStateSpace")
    tf <- tempfile(fileext = ".stan")
    cat("functions{\n#include ssm.stan\n}\nmodel{\n}\n",
        file = tf)
    mc <- rstan::stanc_builder(tf, isystem = isystem)
    expect_is(mc, "list")
    expect_named(mc,
                 c("status", "model_cppname", "cppcode",
                   "model_name", "model_code"))
    capture_output({
      mod <- rstan::stan_model(model_code = mc[["model_code"]])
    })
    expect_is(mod, "stanmodel")
  })

  for (stanfile in dir(system.file("stan/models", package = "StanStateSpace"),
                       full.names = TRUE)) {

    test_that(sprintf("%s parses and compiles", stanfile), {
      isystem <- system.file("stan/include", package = "StanStateSpace")
      mc <- rstan::stanc_builder(stanfile, isystem = isystem)
      expect_is(mc, "list")
      expect_named(mc,
                   c("status", "model_cppname", "cppcode",
                     "model_name", "model_code"))
      capture.output({
        mod <- rstan::stan_model(model_code = mc[["model_code"]])
      })
      expect_is(mod, "stanmodel")
    })
  }
}
