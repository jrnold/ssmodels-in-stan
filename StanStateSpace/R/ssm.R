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

extract_ssm_parameters <- function(extractor, params) {
  f <- function(el, x) {
    d <- dim(x)[1:2]
    aperm(array(aperm(x[ , , el[["start"]]:el[["end"]], drop = FALSE]),
                c(el[["dim"]], rev(d))))
  }
  if (!is.null(params)) {
    extractor <- extractor[params]
  }
  map(extractor, ~ f(.x, x))
}


#' Extract parameters from Stan State Space Model results
#'
#' Extract results from the output of the \code{ssm_filter} function
#' in stan.
#'
#' @param x Array containing results from \code{ssm_filter} Stan function. This
#'   is usually extracted from a \code{stanfit} object.
#' @param m Number of states
#' @param p Dimension of observation vector
#' @param params Parameters to extract from the filter. One of
#'   \code{"v"}, \code{"Finv"}, \code{"K"}, \code{"a"}, \code{"P"}.
#'   If \code{NULL}, all values are extracted.
#' @return A named \code{list} of \code{array} objects, one for each parameter.
extract_from_ssm_filter <- function(x, m = NULL, p = NULL, q = NULL, params = NULL) {
  extractor <- gen_ssm_filter_extractor(m, p)
}
