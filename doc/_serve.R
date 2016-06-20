suppressPackageStartupMessages({
  library("optparse")
  library("bookdown")
})
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--preview"), action = "store_true",
                     default = TRUE,
                     help = "Whether to render the modified/added chapters only [default]")
parser <- add_option(parser, c("-a", "--all"), action = "store_true",
                     dest = "preview",
                     help = "Whether to render the whole book")
pargs <- parse_args(parser, args = commandArgs(TRUE))
print(pargs)

serve_book(dir = ".",
           preview = pargs$preview,
           daemon = TRUE)
