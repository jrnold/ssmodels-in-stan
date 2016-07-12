#' Stan Paths and Files
#'
#' \code{ssm_stan_include_path} returns the path to
#' the directory with the Stan user-defined functions provided
#' by this package. This is used with the \code{isystem}
#' argument in \code{\link[rstan]{stanc_builder}}.
#'
#' \code{ssm_stan_models_path} returns the path to
#' the directory containing the models provided by this package.
#'
#' \code{ssm_available_models} lists the names of the stan models
#' provided by this package.
#'
#' \code{ssm_stan_model} runs a Stan model included in this package.
#' This function is a wrapper for \code{\link[rstan]{stan}}, that
#' sets the correct path to the model and includes the functions provided
#' in this package.
#'
#' @export
ssm_stan_include_path <- function() {
  system.file("stan/include", package = "StanStateSpace")
}

#' @rdname ssm_stan_include_path
#' @export
ssm_stan_model_path <- function() {
  system.file("stan/models", package = "StanStaceSpace")
}

#' @rdname ssm_stan_include_path
#' @export
ssm_available_models <- function() {
  dir(ssm_stan_model_path(), pattern = "\\.stan$")
}

#' @param file The name of the Stan model to run
#' @param ... Passed to \code{\link[rstan]{stan}}
#' @rdname ssm_stan_include_path
#' @export
ssm_stan_model <- function(file, ...) {
  if (!exists(file)) {
    file <- file.path(ssm_stan_model_path(), file)
  }
  stan(model_code = stanc_builder(file, isystem = ssm_stan_include_path()),
       model_name = gsub("\\.stan$", "", basename(file)), ...)
}

gen_ssm_extractor <- function(...) {
  params <- list(...)
  ret <- list()
  start <- 0
  for (i in seq_along(params)) {
    x <- params[[i]]
    if (x[["type"]] == "symmetric_matrix") {
      n <- x[["dim"]][1]
      len <- n * (n + 1) / 2
    } else {
      len <- prod(x[["dim"]])
    }
    ret[[names(params)[i]]] <-
      list(start = start + 1,
           end = start + len,
           len = len,
           dim = x[["dim"]],
           type = x[["type"]])
    start <- start + len
  }
  attr(ret, "vector_length") <- ret[[length(params)]][["end"]]
  ret
}

ssm_extractors <- within(list(), {
  filter <- function(m, p, q) {
    gen_ssm_extractor(loglik = list(dim = 1, type = "real"),
                      v = list(dim = p, type = "vector"),
                      Finv = list(dim = c(p, p), type = "symmetric_matrix"),
                      K = list(dim = c(m, p), type = "matrix"),
                      a = list(dim = m, type = "vector"),
                      P = list(dim = c(m, m), type = "symmetric_matrix"))
  }

  filter_states <- function(m, p, q) {
    gen_ssm_extractor(a = list(dim = m, type = "vector"),
                      P = list(dim = c(m, m), type = "symmetric_matrix"))
  }

  smooth_state <- function(m, p, q) {
    gen_ssm_extractor(alpha = list(dim = m, type = "vector"),
                      V = list(dim = c(m, m), type = "symmetric_matrix"))
  }

  smooth_eps <- function(m, p, q) {
    gen_ssm_extractor(mean = list(dim = p, type = "vector"),
                      var = list(dim = c(p, p), type = "symmetric_matrix"))
  }

  smooth_eta <- function(m, p, q) {
    gen_ssm_extractor(mean = list(dim = q, type = "vector"),
                      var = list(dim = c(q, q), type = "symmetric_matrix"))
  }

  sim_rng <- function(m, p, q) {
    gen_ssm_extractor(y = list(dim = p, type = "vector"),
                      a = list(dim = m, type = "vector"),
                      eps = list(dim = p, type = "vector"),
                      eta = list(dim = q, type = "vector"))
  }
})

ssm_extract_param <- function(param, x) {
  # The dimensions of the array are (iteration, time, returndim)
  if (param[["type"]] == "symmetric_matrix") {
    iter <- dim(x)[1]
    time <- dim(x)[2]
    n <- floor(sqrt(2 * param[["len"]]))
    k <- param[["start"]]:param[["end"]]
    ret <- array(NA_real_, c(iter, time, n, n))
    for (i in seq_len(iter)) {
      for (j in seq_len(time)) {
        ret[i, j, , ] <- vector_to_symmat(as.numeric(x[i, j, k]))
      }
    }
    ret
  } else {
    if (length(param[["dim"]]) == 1) {
      x[ , , param[["start"]]:param[["end"]], drop = FALSE]
    } else {# length == 2
      iter <- dim(x)[1]
      time <- dim(x)[2]
      veclen <- dim(x)[3]
      m <- param[["dim"]][1]
      n <- param[["dim"]][2]
      k <- param[["start"]]:param[["end"]]
      ret <- array(NA_real_, c(iter, time, m, n))
      for (i in seq_len(iter)) {
        for (j in seq_len(time)) {
          ret[i, j, , ] <- matrix(x[i, j, k], m, n)
        }
      }
      ret
    }
  }
}

#' Extract parameters from Stan State Space Model results
#'
#' Extract results from the output of the \code{ssm_filter} function
#' in stan.
#'
#' @param x An array containing results from one of the "ssm" Stan functions:
#'   \code{ssm_filter}, \code{ssm_filter_states}, \code{ssm_smooth_eta},
#'   \code{ssm_smooth_eps}, \code{ssm_smooth_state}, or \code{ssm_sim_rng}.
#'   This array is usually extracted from a \code{stanfit} object using
#'   \code{\link[rstan]{extract}}.
#' @param m Dimension of the state vector
#' @param p Dimension of observation vector
#' @param q Dimension of the state disturbance vector
#' @param type Character. The function which generated the vector.
#' @param params Character. Then names of arameters to extract from the filter. See @details
#'   for the available parameters for each type.
#'   If \code{NULL}, then all parameters for that type are extracted.
#' @return A named \code{list} of \code{array} objects, one for each parameter.
#'
#' @details
#'
#' The parameters returned by each \code{type} are:
#' \describe{
#' \item{\code{"filter"}}{\code{"v"}, \code{"Finv"}, \code{"K"}, \code{"a"}, \code{"P"}}
#' \item{\code{"filter_states"}}{\code{"a"}, \code{"P"}}
#' \item{\code{"smooth_state"}}{\code{"alpha"}, \code{"V"}}
#' \item{\code{"smooth_eta"}}{\code{"mean"}, \code{"var"}}
#' \item{\code{"smooth_eps"}}{\code{"mean"}, \code{"var"}}
#' }
#'
#' @export
ssm_extract <- function(x, m, p, q = m,
                        type = c("filter", "filter_states",
                                 "smooth_state", "smooth_eps",
                                 "smooth_eta", "sim_rng"),
                        params = NULL) {
  # states must be integers >= 0
  assert_that(is.count(m))
  assert_that(is.count(p))
  assert_that(is.count(q))
  # the state innovations must be <= the number of states
  assert_that(q <= m)
  # x should be a numeric array with length 3
  assert_that(is.numeric(x))
  assert_that(length(dim(x)) == 3L)
  # Check that type is one of the supported types
  # It would be better if this was directly tied to ssm_extractor names
  type <- match.arg(type)
  extractor <- ssm_extractors[[type]](m, p, q)
  # Check that vector length is the one expected by the
  # extractor
  nx = attr(extractor, "vector_length")
  if (dim(x)[3] != nx) {
    stop(sprintf("For m = %d, p = %d, q = %d and type = %s, expected dim(x)[3] == %d",
                 m, p, q, type, nx))
  }
  bad_params <- setdiff(params, names(extractor))
  if (length(bad_params) > 0) {
    stop(sQuote("param"), " value(s) ",
         paste(dQuote(bad_params), collapse = ", "),
         " are invalid.",
         sQuote("param"), " for type ", dQuote(type),
         " must be one of ",
         paste(dQuote(names(extractor)), collapse = ", "))
  }
  if (!is.null(params)) {
    extractor <- extractor[params]
  }
  map(extractor, ssm_extract_param, x = x)
}
