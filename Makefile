RSCRIPT = Rscript

pdffile = dlmstan.pdf
srcfile = $(pdffile:%.pdf=%.Rnw)
auxfile = $(pdffile:%.pdf=%.bcf)

all: $(pdffile)

%.pdf: %.tex
	xelatex -interaction nonstopmode $<
	biber $(auxfile)
	xelatex -interaction nonstopmode $<
	xelatex -interaction nonstopmode $<

%.tex: %.Rnw
	$(RSCRIPT) -e 'library(knitr);knit("$<",output="$@")'

