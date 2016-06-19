#' Download data from the Commandeur and Koopman book "An Introduciton to State Space Time Series Analysis"
#'
#'
ZIP_FILE_URL <- "http://staff.feweb.vu.nl/koopman/projects/ckbook/OxCodeAll.zip"
DATAFILES <- paste("OxCodeIntroStateSpaceBook",
                   c(
                     "Chapter_8/NorwayFinland.txt",
                     "Chapter_7/UKinflation.txt"
                   ),
                   sep = "/")
DST <- "data-raw"


#' Unique datasets
# NorwayFinland.txt
# UKdriversKSI.txt
# UKfrontrearseatKSI.txt
# UKinflation.txt
# logUKpetrolprice.txt
# logpetrol.txt
# norway.txt
# signalchapter7-3.txt

dir.create(DST, showWarnings = FALSE)
resp <- httr::GET(ZIP_FILE_URL)
tf <- tempfile()
writeBin(resp$content, tf)
unzip(tf, files = DATAFILES, exdir = DST)

