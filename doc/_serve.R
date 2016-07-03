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
parser <- add_option(parser, c("-D", "--dir"),
                     action = "store",
                     default =  formals(bookdown::serve_book)[["dir"]],
                     dest = "dir",
                     help = "Directory with the sources")
parser <- add_option(parser, c("-o", "--output-dir"),
                     action = "store",
                     type = "character",
                     dest = "output_dir",
                     default = formals(bookdown::serve_book)[["output_dir"]],
                     help = "Output directory")
pargs <- parse_args(parser, args = commandArgs(TRUE))
print(pargs)

serve_book(dir = pargs$dir,
           output_dir = pargs$output_dir,
           preview = pargs$preview,
           in_session = FALSE,
           daemon = TRUE)
