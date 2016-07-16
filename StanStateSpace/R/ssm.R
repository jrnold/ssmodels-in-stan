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
  ret <- vector(length = length(params), mode = "list")
  names(ret) <- names(params)
  # This needs to be done in a loop and sequentially in order to
  # keep track of the start index for each parameter
  start <- 0
  for (i in seq_along(params)) {
    x <- params[[i]]
    par_name <- names(params)[i]
    if (x[["type"]] == "symmetric_matrix") {
      n <- x[["dim"]][1]
      len <- n * (n + 1) / 2
      parnames <- character(len)
      idx_row <- integer(len)
      idx_col <- integer(len)
      .k <- 1
      for (.j in 1:n) {
        for (.i in .j:n) {
          parnames[.k] <- sprintf("%s[%d,%d]", par_name, .i, .j)
          idx_row[.k] <- .i
          idx_col[.k] <- .j
          .k <- .k + 1
        }
      }
    } else if (x[["type"]] == "vector") {
      len <- x[["dim"]][1]
      parnames <- sprintf("%s[%s]", par_name, seq_len(len))
      idx_row <- seq_len(len)
      idx_col <- rep(1L, len)
    } else if (x[["type"]] == "matrix") {
      .rows <- x[["dim"]][1]
      .cols <- x[["dim"]][2]
      len <- .rows * .cols
      parnames <- paste0(par_name, "[",
                         apply(expand.grid(seq_len(.rows), seq_len(.cols)), 1,
                               paste0, collapse = ","), "]")
      idx_row <- rep(seq_len(.rows), .cols)
      idx_col <- rep(seq_len(.cols), each = .rows)
    } else if (x[["type"]] == "real") {
      len <- 1L
      parnames <- names(params)[i]
      idx_row <- NA_integer_
      idx_col <- NA_integer_
    } else {
      stop("type = ", dQuote(x[["type"]]), " not recognized.")
    }
    ret[[par_name]] <-
      list(start = start + 1,
           end = start + len,
           idx = start + seq_len(len),
           len = len,
           dim = x[["dim"]],
           type = x[["type"]],
           parindex = cbind(idx_row, idx_col),
           parnames = parnames)
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
    k <- param[["idx"]]
    ret <- array(NA_real_, c(iter, time, n, n))
    for (i in seq_len(iter)) {
      for (j in seq_len(time)) {
        ret[i, j, , ] <- vector_to_symmat(as.numeric(x[i, j, k]))
      }
    }
    ret
  } else {
    if (length(param[["dim"]]) == 1) {
      x[ , , param[["idx"]], drop = FALSE]
    } else {# length == 2
      iter <- dim(x)[1]
      time <- dim(x)[2]
      m <- param[["dim"]][1]
      n <- param[["dim"]][2]
      k <- param[["idx"]]
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
    stop(sprintf(paste0("For m = %d, p = %d, q = %d and type = %s,",
                        "expected dim(x)[3] == %d"),
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

#' Extract State Space Parameters from a Stanfit summary
#'
#' @param x Result of the \code{\link[rstan:summary,stanfit-method]{summary}} method for a \code{\link[rstan]{stanfit}} object.
#' @param par Name of the parameter containing state space parameters
#' @param m Number of states
#' @param p Size of observation vector in the state space model
#' @param q Size of the state disturbance vector in the state space model
#' @param type The name of the state space function which generated the vector.
#' @param chains If \code{FALSE}, then use the summary over all chains. If \code{TRUE}, then use the individual chain summaries.
#' @return A \code{data_frame} with parameters as rows, and summary statistics and metadata about the parameters in the columns:
#'    \describe{
#'      \item{\code{par_id}}{Parameter identifier, e.g. \code{"Finv[1,2]"}}
#'      \item{\code{parameter}}{Parameter name, e.g. \code{"Finv"}}
#'      \item{\code{index}}{Index number in the original vector.}
#'      \item{\code{index_row}}{Row index of this element within the parameter}
#'      \item{\code{index_col}}{Column index of this element within the parameter}
#'      \item{\code{mean}}{Mean}
#'      \item{\code{se_mean}}{Standard error of the mean.  Only if \code{chains = FALSE}}
#'      \item{\code{sd}}{Standard deviation}
#'      \item{\code{10\%}}{10th percentile}
#'      \item{\code{90\%}}{90th percentile}
#'      \item{\code{n_eff}}{Number of effective samples.  Only if \code{chains = FALSE}}
#'      \item{\code{Rhat}}{R-hat. Only if \code{chains = FALSE}}
#'    }
#'
#' @export
ssm_extract_summary <- function(x, par, m, p, q = m,
                                type = c("filter", "filter_states",
                                         "smooth_state", "smooth_eps",
                                         "smooth_eta", "sim_rng"),
                                chains = FALSE) {
  # states must be integers >= 0
  one_of <- NULL
  assert_that(is.count(m))
  assert_that(is.count(p))
  assert_that(is.count(q))
  # the state innovations must be <= the number of states
  assert_that(q <= m)
  # Check that type is one of the supported types
  # It would be better if this was directly tied to ssm_extractor names
  type <- match.arg(type)
  pattern <- sprintf("^%s\\[(\\d+),(\\d+)\\]$", par)
  extractor <- ssm_extractors[[type]](m, p, q)
  parameters <- map_df(extractor, function(.) {
    data_frame(parameter = .[["parnames"]],
               index = .[["idx"]],
               index_row = .[["parindex"]][ , 1],
               index_col = .[["parindex"]][ , 2])
  }, .id = "par_id")
  if (!chains) {
    param_rows <- str_detect(rownames(x[["summary"]]), pattern)
    parnames <- rownames(x[["summary"]])[param_rows]
    dat <- as_data_frame(x[["summary"]][param_rows, ])
    dat[["par_id"]] <- parnames
    dat <- separate_(dat, "par_id", c(".orig", "time", "index"),
                    remove = FALSE, extra = "drop", convert = TRUE)
    dat <- select(dat, -one_of(c(".orig", "par_id")))
  } else {
    parnames <- dimnames(x[["c_summary"]])[1]
    variables <- dimnames(x[["c_summary"]])[2]
    param_rows <- str_detect(parnames, pattern)
    f <- function(.x, i, parnames, variables) {
      .df <- as_data_frame(.x[i, ])
      colnames(.df) <- variables
      .df[["par_id"]] <- parnames
      .df <- separate_(.df, "par_id", c(".orig", "time", "index"),
                      remove = FALSE, extra = "drop", convert = TRUE)
      .df <- select(.df, -one_of(c(".orig", "par_id")))
      .df
    }
    dat <- map_df(array_branch(x[["c_summary"]], 3),
           f, i = param_rows, parnames = parnames, variables = variables,
           .id = "chain")
  }
  nx = attr(extractor, "vector_length")
  if (max(dat[["index"]]) != nx) {
    stop(sprintf(paste("For m = %d, p = %d, q = %d and type = %s,",
                       "expected the number of parameters to equal %d"),
                 m, p, q, type, nx))
  }
  left_join(parameters, dat, by = "index")
}

#' Tidy the results of Stanfit summary objects
#'
#' @param x The list returned by calling \code{summary} on a \code{Stanfit} object.
#' @return A list with two elements,
#' \describe{
#' \item{all}{A data frame with summary statistics for all the chains combined.}
#' \item{chains}{A data frame with summary statistics for each chain.}
#' }
#'
#' @export
tidy_stan_summary <- function(x) {
  params <- parse_stan_parnames(rownames(x[["summary"]]))
  ret <- list()
  ret[["all"]] <- bind_cols(params, as_data_frame(x[["summary"]]))
  ret[["chains"]] <-
    map_df(array_branch(x[["c_summary"]], 3),
           function(.data) bind_cols(params, as_data_frame(.data)),
           .id = "chain")
  ret
}