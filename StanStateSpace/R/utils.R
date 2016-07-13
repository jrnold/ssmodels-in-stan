#' @importFrom purrr map map_df array_branch
#' @importFrom stringr str_detect
#' @importFrom tibble rownames_to_column as_data_frame
#' @importFrom dplyr data_frame left_join
#'
#' @import stats
#'
#' @import assertthat
#' @import lubridate
#' @import rstan
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
#' @param rows logical. If \code{TRUE}, then the selection matrix
#'   selects rows, if false, then it selects columns.
#'
#' @param k Index of the diagonal. 0 is the main diagonal, positive numbers for
#'    upper diagonals, negative numbers for lower diagonals.
#' @return For \code{unit_lower_tri} and \code{unit_upper_tri} a
#'   matrix of dimensions \code{c(m, n)}. For \code{selection_matrix}
#'   a matrix with dimensions \code{m, length(idx)} if \code{rows = TRUE},
#'   else a matrix with dimensions \code{c(length(idx), m}
#'
#' @export
unit_lower_tri <- function(m, n = m, diag = TRUE) {
  matrix(as.numeric(lower.tri(matrix(0, m, n), diag = diag)), m, n)
}


#' @describeIn unit_lower_tri Generate an unit upper triangular matrix.
#' @export
unit_upper_tri <- function(m, n = m, diag = TRUE){
  matrix(as.numeric(upper.tri(matrix(0, m, n), diag = diag)), m, n)
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


#' @describeIn unit_lower_tri Generate matrix with ones on the diagonal
#' @export
eye <- function(m, n = m, k = 0) {
  if (k == 0) {
    x <- diag(1, m, n)
  } else {
    # not identity
    x <- matrix(0, m, n)
    if (k < 0 && abs(m > abs(k))) {
      for (i in 1:(m + k)) {
        x[i - k, i] <- 1
      }
    } else if (k > 0 && abs(n > k)) {
      for (i in 1:(n - k)) {
        x[i, i + k] <- 1
      }
    }
  }
  x
}


#' Vector to symmetric matrix
#'
#' @param x vector with n (n + 1) / 2 elements.
#' @return A n x n symmetric matrix
#'
#' This fills in the matrix assuming that the elements
#' are from the lower triangular part of the symmetric
#' matrix, and were filled in column-wise.
#'
#' @export
vector_to_symmat <- function(x) {
  # This should be rewritten in RCpp. It is not very
  # idiomatic at the moment.
  num_el <- length(x)
  n <- floor(sqrt(2 * num_el))
  # use x[1] to ensure newmat has the same type as x
  newmat <- matrix(x[1], n, n)
  k <- 1
  for (j in 1:n) {
    for (i in j:n) {
      newmat[i, j] <- x[k]
      if (i != j) {
        newmat[j, i] <- x[k]
      }
      k <- k + 1
    }
  }
  newmat
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
    daze <- floor(tx %% 1 * (365 + leap_year(yr)))
    dates <- yr + days(daze)
  }
  dates
}
