library("httr")
library("zoo")
library("dplyr")
library("tibble")
library("tidyr")

#' Tussell (2011)
#'
#' Download Tussell data divsas.rda and save as a data frame
#'
get_divisas <- function() {
  dest <- "data/divisas.rda"
  if (!file.exists(dest)) {
    con <- url("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v039i02/divisas.rda")
    load(con)
    close(con)
    colnames(divisas) <- c("AUD","BEF","CHF","DEM",
                           "DKK","ESP","FRF","GBP","ITL",
                           "JPY","NGL","SEK")
    divisas <- as.data.frame(divisas) %>%
      rownames_to_column(var = "date") %>%
      mutate(date = as.Date(date, "%Y-%m-%d")) %>%
      gather(currency, fxrate, -date)
    save(divisas, file = dest)
    message("Created ", dest)
  } else {
    message(dest, " already exists. Nothing done.")
  }
}


build <- function(clean = FALSE) {
  if (clean) {
    message("Deleting files in data/")
    for (filename in dir("data", all.files = FALSE, full.names = TRUE)) {
      file.remove(filename)
      message("Deleted ", filename)
    }
  }
  get_divisas()
}