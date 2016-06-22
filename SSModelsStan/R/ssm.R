#' @importFrom purrr map array_tree
NULL

gen_ssm_extractor <- function(...) {
  params <- list(...)
  ret <- list()
  start <- 0
  for (i in seq_along(params)) {
    x <- params[[i]]
    len <- prod(x)
    ret[[names(params)[i]]] <-
      list(start = start + 1,
           end = start + len,
           len = len,
           dim = x)
    start <- start + len
  }
  ret
}

gen_ssm_filter_extractor <- function(m, p) {
  gen_ssm_extractor(loglik = 1,
                    v = p,
                    Finv = c(p, p),
                    K = c(m, p),
                    a = m,
                    P = c(m, m))
}

#' Extract parameters from Stan State Space results
#'
#'
#' Extract results from the output of the \code{ssm_filter} function
#' in stan.
#'
#'
#' @param x Array containing results from \code{ssm_filter} Stan function. This
#'   is usually extracted from a \code{stanfit} object.
#' @param m Number of states
#' @param p Dimension of observation vector
#' @param params Parameters to extract from the filter. One of
#'   \code{"v"}, \code{"Finv"}, \code{"K"}, \code{"a"}, \code{"P"}.
#'   If \code{NULL}, all values are extracted.
#' @return A named \code{list} of \code{array} objects, one for each parameter.
#'
extract_param_from_ssm_filter <- function(x, m, p, params = NULL) {
  f <- function(x, extractor) {
    g <- function(x, el) {
      .n <- nrow(x)
      array(t(x[ , el[["start"]]:el[["end"]]]), c(el[["dim"]], .n))
    }
    map(extractor, ~ g(x, .x))
  }
  extractor <- gen_ssm_filter_extractor(m, p)
  if (!is.null(params)) {
    extractor <- extractor[params]
  }
  map(array_tree(x, 1), ~ f(.x, extractor))
}
