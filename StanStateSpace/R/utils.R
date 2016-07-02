#' @importFrom purrr map array_tree
#' @import stats
#' @import lubridate
NULL

#' Upper and lower unit triangle matrices
#'
#' Some useful functions for generating matrices commonly encountered
#' in state space models.
#'
#'
#' @param m integer. Number of rows in the matrix
#' @param n integer. Number of columns in the matrix
#' @param diag logical. Include diagonal?
#' @param idx integer. Indexes of rows (if \code{rows = TRUE}),
#'   or columns (if \code{rows = FALSE}) to select.
#' @param j The index to set to one.
#' @param rows logical. If \code{TRUE}, then the selection matrix
#'   selects rows, if false, then it selects columns. In the
#' @return For \code{unit_lower_tri} and \code{unit_upper_tri} a
#'   matrix of dimensions \code{c(m, n)}. For \code{selection_matrix}
#'   a matrix with dimensions \code{m, length(idx)} if \code{rows = TRUE},
#'   else a matrix with dimensions \code{c(length(idx), m}
#'
#' @export
unit_lower_tri <- function(m, n = m, diag = TRUE) {
  as.numeric(lower.tri(matrix(0, m, n), diag = diag))
}

#' @describeIn unit_lower_tri Generate an unit upper triangular matrix.
#' @export
unit_upper_tri <- function(m, n = m, diag = TRUE){
  as.numeric(upper.tri(matrix(0, m, n), diag = diag))
}

#' @describeIn unit_lower_tri Generate a selection matrix.
#' @export
selection_matrix <- function(m, idx, rows = TRUE) {
  n <- length(idx)
  x <- matrix(0, n, m)
  for (i in seq_along(idx)) {
    x[i, idx[i]] <- 1
  }
  if (!rows) {
    x <- t(x)
  }
  x
}
#' @describeIn unit_lower_tri Generate selection vector
#' @export
eye <- function(n, j = 1) {
  x <- rep(0, n)
  x[j] <- 1
  x
}

#' Convert
#'
#' This function assumes that the data in the \code{ts} object are
#' sampled daily. There case when \code{frequency = 12} is treated as
#' monthly frequency, and \code{frequency = 4} is treated as
#' quarterly frequency.
#'
#' @param x A \code{ts} object. The \code{tsp} attribute is assumed to have
#'   be sampled at the daily frequency with integer parts representing years.
#' @param ... Not used
#' @return A \code{Date} vector
#'
#' @export
as.Date.ts <- function(x, ...) {
  f <- frequency(x)
  tx <- as.numeric(time(x))
  yr <- ymd(paste(floor(tx), "01", "01"))
  if (f == 12) {
    dates <- yr %m+% months(round((tx %% 1) * 12))
  } else if (f == 4) {
    dates <- yr %m+% months(round((tx %% 1) * 4 * 3))
  } else if (f == 1) {
    dates <- yr
  } else {
    daze <- tx %% 1 * (365 + leap_year(yr))
    dates <- yr + days(daze)
  }
  dates
}
