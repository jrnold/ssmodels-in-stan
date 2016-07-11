library("httr")
library("zoo")
library("dplyr")
library("tibble")
library("tidyr")
library("readr")
library("stringr")

run_if_not_exists <- function(obj, expr) {
  envir <- parent.frame()
  dest <- file.path("data", paste0(obj, ".rda"))
  if (!file.exists(dest)) {
    assign(obj, eval(expr, envir = envir))
    save(obj, file = dest)
    message("Created ", dest)
  } else {
    message(dest, " already exists. Nothing done.")
  }
}

#' Tussell (2011)
#'
#' Download Tussell data divsas.rda and save as a data frame
#'
get_divisas <- function() {
  run_if_not_exists("divisas", {
    con <- url("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v039i02/divisas.rda")
    load(con)
    close(con)
    colnames(divisas) <- c("AUD","BEF","CHF","DEM",
                           "DKK","ESP","FRF","GBP","ITL",
                           "JPY","NGL","SEK")
    as.data.frame(divisas) %>%
      rownames_to_column(var = "date") %>%
      mutate(date = as.Date(date, "%Y-%m-%d")) %>%
      gather(currency, fxrate, -date)
  })
}

get_capm <- function() {
  run_if_not_exists("capm", {
    tf <- tempfile()
    r <- GET("https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v041i04/P.txt",
             write_disk(tf))
    capm <- read.table(tf) %>%
      rownames_to_column("yq") %>%
      separate(yq, c("year", "month")) %>%
      mutate(date = as.Date(paste(year, "-", month, "-01", sep = ""), "%Y-%m-%d")) %>%
      select(-year, -month) %>%
      select(date, everything()) %>%
      arrange(date)
    names(capm) <- str_to_lower(names(capm))
    capm
  })
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
  get_capm()
}