library("rstan")
devtools::install("../../StanStateSpace")

# Run a stan test
# @param file name of function. The file is "stan/test_{file}.stan"
# @param init Inputs to the function
# @param data Other data needed to initialize the model
# @return An array with the output of the function.
test_stan_function <- function(FUN, data = NULL, init = NULL,
                               output = TRUE) {
  filename <- file.path( "build",  paste0("test_", FUN))
  system2("make", args = c(filename))
  out_tmpfile <- tempfile(fileext = ".csv")
  args <- c("sample",
            "num_samples=1",
            "num_warmup=0",
            "algorithm=fixed_param",
            "output",
            paste0("file=", out_tmpfile))
  if (!is.null(data)) {
    data_tmpfile <- tempfile(fileext = ".R")
    rstan::stan_rdump(names(data), file = data_tmpfile,
                      envir = as.environment(data))
    args <- append(args, c("data", paste0("file=", data_tmpfile)))
  }
  if (!is.null(init)) {
    init_tmpfile <- tempfile(fileext = ".R")
    rstan::stan_rdump(names(init), file = data_tmpfile,
                      envir = as.environment(data))
    args <- append(args, paste0("init=", init_tmpfile))
  }
  msg <- capture.output({
    rc <- devtools::system_check(filename, args = args)
  })
  # If error, then raise it
  if (rc > 0) {
    stop(msg)
  }
  if (output) {
    # IF output, then return data
    rstan::read_stan_csv(out_tmpfile)
  } else {
    # else return the return code
  }
}
