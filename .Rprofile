suppressPackageStartupMessages({
  library("ggplot2")
  library("rstan")
  library("dplyr")
  library("lubridate")
  library("rprojroot")
})
devtools::install("SSModelsStan")

# Convert ts object to dates assuming that it has start and end years
ts_to_date <- function(x) {
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
