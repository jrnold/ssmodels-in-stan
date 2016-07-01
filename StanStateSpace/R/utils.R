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

# setClass("SSM",
#          representation(m = "integer", p = "integer", r = "integer"))
#
# ssm_poly <- function(order = 1, H, R, Q, a1, P1) {
#   m <- order
#   p <- 1
#   r <- ncol(R)
#   ret <- list()
#   ret[["Z"]] <- diag(1, p, m)
#   ret[["H"]] <- H
#   ret[["T"]] <- unit_upper_tri(m, m)
#   ret[["R"]] <- R
#   ret[["Q"]] <- Q
#   ret[["a1"]] <- a1
#   ret[["P1"]] <- P1
#   ret[["c"]] <- rep(0, m)
#   ret[["d"]] <- rep(0, p)
#   ret[["m"]] <- m
#   ret[["p"]] <- p
#   ret[["r"]] <- r
#   class(ret) <- c("ssm_constant", "ssm")
# }
#
# ssm_seasonal <- function(frequency, sigma_eta, R, Q, a1, P1) {
#   m <- frequency - 1L
#   p <- 1
#   r <- ncol(R)
#   ret <- list()
#   ret[["Z"]] <- diag(1, p, m)
#   ret[["H"]] <- H
#   TT <- matrix(0, m, m)
#   TT[1, ] <- -1
#   for (i in 2:m) {
#     TT[i, i - 1] <- 1
#   }
#   ret[["R"]] <- diag(1, m, 1)
#   ret[["Q"]] <- matrix(sigma_eta^2, 1, 1)
#   ret[["a1"]] <- a1
#   ret[["P1"]] <- P1
#   ret[["m"]] <- m
#   ret[["p"]] <- p
#   ret[["r"]] <- 1
#   class(ret) <- c("ssm_constant", "ssm")
# }
#
