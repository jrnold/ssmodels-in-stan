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
  ret
}

gen_ssm_filter_extractor <- function(m, p) {
  gen_ssm_extractor(loglik = list(dim = 1, type = "real"),
                    v = list(dim = 1, type = "vector"),
                    Finv = list(dim = c(p, p), type = "symmetric_matrix"),
                    K = list(dim = c(m, p), type = "matrix"),
                    a = list(dim = m, type = "vector"),
                    P = list(dim = c(m, m), type = "symmetric_matrix"))
}

gen_ssm_filter_states_extractor <- function(m) {
  gen_ssm_extractor(a = list(dim = m, type = "vector"),
                    P = list(dim = c(m, m), type = "symmetric_matrix"))
}

gen_ssm_smooth_state_extractor <- function(m) {
  gen_ssm_extractor(alpha = list(dim = m, type = "vector"),
                    V = list(dim = c(m, m), type = "symmetric_matrix"))
}

gen_ssm_smooth_eps_extractor <- function(p) {
  gen_ssm_extractor(mean = list(dim = p, type = "vector"),
                    var = list(dim = c(p, p), type = "symmetric_matrix"))
}

gen_ssm_smooth_eta_extractor <- function(q) {
  gen_ssm_extractor(mean = list(dim = q, type = "vector"),
                    var = list(dim = c(q, q), type = "symmetric_matrix"))
}

gen_ssm_sim_rng_extractor <- function(m, p, q) {
  gen_ssm_extractor(y = list(dim = p, type = "vector"),
                    a = list(dim = m, type = "vector"),
                    eps = list(dim = p, type = "vector"),
                    eta = list(dim = q, type = "vector"))
}

extract_ssm_param <- function(param, x) {
  # The dimensions of the array are (iteration, time, returndim)
  if (param[["type"]] == "symmetric_matrix") {
      x[ , , param[["start"]]:param[["end"]], drop = FALSE],
  } else {
    if (length(param[["dim"]]) == 1) {
      aperm(x[ , , param[["start"]]:param[["end"]], drop = FALSE])
    } else {# length == 2
      d <- dim(x)[1:2]
      aperm(array(aperm(x[ , , param[["start"]]:param[["end"]], drop = FALSE]),
                  c(param[["dim"]], rev(d))))
    }
  }
}

extract_ssm_all_params <- function(extractor, params) {
  if (!is.null(params)) {
    extractor <- extractor[params]
  }
  #map(extractor, ~ f(.x, x))
}


#' Extract parameters from Stan State Space Model results
#'
#' Extract results from the output of the \code{ssm_filter} function
#' in stan.
#'
#' @param x Array containing results from \code{ssm_filter} Stan function. This
#'   is usually extracted from a \code{stanfit} object.
#' @param m Dimension of the state vector
#' @param p Dimension of observation vector
#' @param q Dimension of the state disturbance vector
#' @param params Parameters to extract from the filter. One of
#'   \code{"v"}, \code{"Finv"}, \code{"K"}, \code{"a"}, \code{"P"}.
#'   If \code{NULL}, all values are extracted.
#' @return A named \code{list} of \code{array} objects, one for each parameter.
# extract_from_ssm_filter <- function(x, m = NULL, p = NULL, q = NULL, params = NULL) {
#   extractor <- gen_ssm_filter_extractor(m, p)
# }
