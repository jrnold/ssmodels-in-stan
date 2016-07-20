library("rstan")
devtools::install("../../StanStateSpace")
library("StanStateSpace")

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

  }
}

#' Tolerance needs to be a little lower due to round tripping.
TOL <- 1e-5

# Positive definite matrix: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.html
rand_pdmat <- function(n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol = n, rnorm(n ^ 2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  t(O) %*% diag(ev) %*% O
}

#' Transform vector of real numbers to stationary AR(p) coefficients
#'
#' Extracted from `stats::arima()`
ar_trans <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI). I just need to worry about one set of coefficients
  # so ignore the other parts
  .Call(stats:::C_ARIMA_transPars, as.numeric(par), as.integer(c(length(par), rep(0, 5))), TRUE)[[1]]
}

#' Transform PACF to ACF
#'
#' Extracted from `stats::arima()`
pacf_to_acf <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI)
  # Undo the tanh transformation in C_ARIMA_tranPars to get back to partial autocorrelation
  ar_trans(atanh(par))
}

#' Transform vector of stationary AR(p) coefficients to real numbers
#'
#' Extracted from `stats::arima()`
ar_invtrans <- function(par) {
  # ARMA object (AR, MA, SAR, SMA, S period, I, SI)
  # goes from AR to pacf, pacf (-1, 1) to (-infty, infty)
  .Call(stats:::C_ARIMA_Invtrans, as.numeric(par), as.integer(c(length(par), rep(0, 5))))
}

#' Transform ACF to PACF
acf_to_pacf <- function(par) {
  # Redo the tanh transformation to get to partial autocorrelations
  tanh(ar_invtrans(par))
}

#' Stationary ARMA covariance
#'
#' Extracted parts from the function `stats::arima`.
arma_init <- function(theta, phi, method = "Gardner1980", tol = 0) {
  if (method == "Gardner1980") {
    Q0 <- .Call(stats:::C_getQ0, as.numeric(phi), as.numeric(theta))
  } else {
    Q0 <- .Call(stats:::C_getQ0bis, as.numeric(phi), as.numeric(theta), tol = tol)
  }
  Q0
}
