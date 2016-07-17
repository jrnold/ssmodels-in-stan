devtools::install("../StanStateSpace")
suppressPackageStartupMessages({
  library("stats")
  library("ggplot2")
  library("rstan")
  library("dplyr")
  library("lubridate")
  library("rprojroot")
  library("StanStateSpace")
})
filter <- dplyr::filter

set.seed(151235)
options(digits = 3)

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  cache = TRUE
)

options(dplyr.print_min = 6, dplyr.print_max = 6)
