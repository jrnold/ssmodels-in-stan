#' Generate indices for all dimensions
#'
#' Create matrix if all indices for a given
#' dimension vector.
#'
#' @param d Array dimensions
#' @return \code{matrix} with dimensions \code{c(prod(dim), dim)}.
#' @keywords internal
expand_grid_dim <- function(d) {
  as.matrix(expand.grid(lapply(as.integer(d), seq_len)))
}

# parnames_from_flatpars <- function(x, pre, sep, post) {
#   FUN <- function(flatpar, parameter, idx, scalar) {
#     mcmc_parnames_pattern_idx(parameter, as.integer(str_split(idx, ",")[[1]]),
#                               scalar, pre, sep, post)
#   }
#   as.character(dlply(x, "flatpar", splat(FUN)))
# }

#' @rdname mcmc_parnames
#' @aliases mcmc_parname_pattern
#' @title Create paramter names for flattened parameters
#'
#' @description Given a parameter name and a matrix of index values, generate
#' names for the unlisted parameters.
#'
#' \describe{
#' \item{\code{mcmc_parnames_bugs}}{Writes BUGS/JAGS style flat parameter names, e.g. \code{"alpha[1,2]"}}.
#' \item{\code{mcmc_parnames_stan}}{Writes Stan style flat parameter names, e.g. \code{"alpha.1.2"}}.
#' \item{\code{mcmc_parnames_underscore}}{Writes parameter names with indexes seperated by underscores, e.g. \code{"alpha_1_2"}}.
#' \item{\code{mcmc_parnames_pattern}}{Writes parameter names with arbitrary patterns.}
#' }
#'
#' @param x \code{character} Parameter name.
#' @param d \code{integer} Dimension of the array.
#' @param pre \code{character} String to put before indices.
#' @param sep \code{character} String used to seperate indices.
#' @param post \code{character} String to put after indices.
#' @param colmajor \code{logical}. If \code{TRUE}, then indices are
#' in column-major order (R's default), else row-major.
#' @return \code{character} vector of flat parameter names.
#'
#' @examples
#' mcmc_parnames_bugs("alpha", c(1, 2))
#' mcmc_parnames_stan("alpha", c(1, 2))
#' mcmc_parnames_underscore("alpha", c(1, 2))
#' mcmc_parnames_pattern("alpha", c(1, 2), "<", ";", ">")
# mcmc_parnames_pattern <- function(x, d, pre=".", sep=".", post="",
#                                   colmajor = TRUE) {
#   scalar <- identical(d, 1L)
#   FUN <- function(i) {
#     mcmc_parnames_pattern_idx(x, i, scalar, pre = pre, sep = sep, post = post,
#                               colmajor = colmajor)
#   }
#   apply(expand_grid_dim(d), 1, FUN)
# }
#
# mcmc_parnames_pattern_idx <- function(x, idx, scalar = FALSE,
#                                       pre = ".", sep = ".",
#                                       post = "", colmajor = TRUE) {
#   if (scalar) {
#     x
#   } else {
#     if (!colmajor) {
#       if (is.matrix(idx)) {
#         idx <- apply(idx, 1, rev)
#       } else {
#         idx <- rev(idx)
#       }
#     }
#     if (is.matrix(idx)) {
#       idxstr <- apply(idx, 1, paste, collapse = sep)
#     } else {
#       idxstr <- paste(idx, collapse = sep)
#     }
#     paste0(x, pre, idxstr, post)
#   }
# }
#
# #' @rdname mcmc_parnames
# #' @aliases mcmc_parnames_stan
# mcmc_parnames_stan <- function(x, d, colmajor = TRUE) {
#   mcmc_parnames_pattern(x, d, ".", ".", "", colmajor)
# }
#
# #' @rdname mcmc_parnames
# #' @aliases mcmc_parnames_bugs
# mcmc_parnames_bugs <- function(x, d, colmajor = TRUE) {
#   mcmc_parnames_pattern(x, d, "[", ",", "]", colmajor)
# }
#
# #' @rdname mcmc_parnames
# #' @aliases mcmc_parnames_underscore
# mcmc_parnames_underscore <- function(x, d, colmajor = TRUE) {
#   mcmc_parnames_pattern(x, d, "_", "_", "", colmajor)
# }
#
# #' Convert between BUGS/JAGS and Stan style flat parameter names
# #'
# #' Utility functions to convert flat parameter names between the BUGS (\code{"alpha[1,1]"}) and
# #' Stan style (\code{"alpha.1.1"}).
# #'
# #' @param x \code{character} vector of flat parameter names.
# #' @return \code{character} vector of the converted flat parameter names.
# #' @export
# #' @examples
# #' stan_parnames <- c("alpha", "beta.1", "gamma.1.1")
# #' bugs_parnames <- c("alpha", "beta[1]", "gamma[1,1]")
# #' identical(bugs_to_stan_parnames(bugs_parnames), stan_parnames)
# #' identical(stan_to_bugs_parnames(stan_parnames), bugs_parnames)
# bugs_to_stan_parnames <- function(x) {
#   parnames_from_flatpars(mcmc_parparser_bugs(x), ".", ".", "")
# }
#
# #' @rdname bugs_to_stan_parnames
# #' @aliases stan_to_bugs_parnamse
# #' @export
# stan_to_bugs_parnames <- function(x) {
#   parnames_from_flatpars(mcmc_parparser_stan(x), "[", ",", "]")
# }
#
# #' @rdname mcmc_parparsers
# #' @title Parse MCMC parameter names
# #'
# #' @description Functions that parse a vector of flat parameter names
# #' and return an object of class \code{\linkS4class{McmcdbFlatpars}}.
# #'
# #' \describe{
# #' \item{\code{mcmc_parparser_stan}}{Parses parameter names
# #' treating each parameter as a scalar. E.g. \code{"beta.1"}
# #' and \code{"beta.2"} will be treated two parameter arrays of
# #' size 1.}
# #' \item{\code{mcmc_parparser_stan}}{Parses parameter names
# #' in the Stan style, e.g. \code{"beta.1.1"}}
# #' \item{\code{mcmc_parparser_guess}}{Tries to guess the format of the parameters}
# #' \item{\code{mcmc_parparser_pattern}}{Parses parameter names using arbitrary patterns.}
# #' }
# #'
# #' @param x \code{character} vector with flat parameter names.
# #' @param pre \code{character} Pattern between parameter name and indices. If a pattern
# #' grouping must be used, use "(?: )".
# #' @param sep \code{character} Pattern seperating each index.
# #' @param post \code{character} Pattern following the indices.
# #' @param colmajor \code{logical}. If \code{TRUE}, then indices are
# #' in column-major order (R's default), else row-major.
# #' @return Object of class \code{McmcdbFlatpars}
# #'
# #' @examples
# #' mcmc_parparser_bugs(c("beta[1,1]", "beta[1,2]"))
# #' mcmc_parparser_stan(c("beta.1.1", "beta.1.2"))
# #' mcmc_parparser_underscore(c("beta_1_1", "beta_1_2"))
# #' mcmc_parparser_pattern(c("beta<1;1>", "beta<1;2>"), "<", ";", ">")
# #' mcmc_parparser_guess(c("beta[1,1]", "beta[1,2]"))
# #' mcmc_parparser_guess(c("beta.1.1", "beta.1.2"))
# #' mcmc_parparser_scalar(c("beta[1,1]", "beta[1,2]"))
# #' # for pattern groups, you must use (?:
# #' mcmc_parparser_pattern(c("beta<1;1>", "beta.1,2"), "[<.]", "[;,]", "(?:>|)")
# mcmc_parparser_scalar <- function(x) {
#   McmcdbFlatpars(data.frame(flatpar = x,
#                             parameter = x,
#                             idx = as.character("1"),
#                             scalar = TRUE,
#                             stringsAsFactors = FALSE))
# }
#
#
# #' @rdname mcmc_parparsers
# #' @aliases mcmc_parparser_pattern
# mcmc_parparser_pattern <- function(x, pre, sep, post, colmajor=TRUE) {
#   regexp <-
#     sprintf("^(.*?)(%s(\\d+(%s\\d+)*)%s)?$", pre, sep, post)
#   x_split <- data.frame(str_match(x, regexp)[ , c(1, 2, 4)],
#                         stringsAsFactors = FALSE)
#   names(x_split) <- c("flatpar", "parameter", "idx")
#   x_split$scalar <- (x_split$idx == "")
#   x_split$idx[x_split$scalar] <- "1"
#   x_split$idx <- gsub(sep, ",", x_split$idx)
#   # If row-major reverse the order of the index
#   if (!colmajor) {
#     x_split$idx <-
#       sapply(str_split(x_split$idx, fixed(",")),
#              function(string) paste(rev(string), collapse=","))
#   }
#   McmcdbFlatpars(x_split)
# }
#
# #' @rdname mcmc_parparsers
# #' @aliases mcmc_parparser_stan
# mcmc_parparser_stan <- function(x) {
#   mcmc_parparser_pattern(x, "\\.", "\\.", "")
# }
#
# #' @rdname mcmc_parparsers
# #' @aliases mcmc_parparser_bugs
# mcmc_parparser_bugs <- function(x) {
#   mcmc_parparser_pattern(x, "\\[", ",", "\\]")
# }
#
# #' @rdname mcmc_parparsers
# #' @aliases mcmc_parparser_underscore
# mcmc_parparser_underscore <- function(x) {
#   mcmc_parparser_pattern(x, "_", "_", "")
# }
#
# has_bracket_index <- function(x) {
#   str_detect(x, "(\\[(\\d+(,\\d+)?)\\])$")
# }
#
# has_dots_index <- function(x) {
#   str_detect(x, "(\\.(\\d+(\\.\\d+)?))$")
# }
#
# #' @rdname mcmc_parparsers
# #' @aliases mcmc_parparser_guess
# mcmc_parparser_guess <- function(x) {
#   if (any(has_bracket_index(x))) {
#     mcmc_parparser_bugs(x)
#   } else if (any(has_dots_index(x))) {
#     mcmc_parparser_stan(x)
#   } else {
#     mcmc_parparser_scalar(x)
#   }
# }