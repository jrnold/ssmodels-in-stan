unit_lower_tri <- function(m, n = m) {
  as.numeric(lower.tri(matrix(0, m, n)))
}
unit_upper_tri <- function(m, n = m){
  as.numeric(upper.tri(matrix(0, m, n)))
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
