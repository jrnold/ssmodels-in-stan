PDF_FILE = ssmodels_in_stan.pdf
HTML_FILE = ssmodels.html
BIB = ../bib/jrnold.bib

build: pdf html

pdf: $(PDF_FILE)

html: $(HTML_FILE)

$(PDF_FILE): $(PDF_FILE:%.pdf=%.tex)
	latexmk -pdf -xelatex $^

$(HTML_FILE): $(HTML_FILE:%.html=%.md)
	pandoc -s --mathjax -f markdown+pipe_tables -o $@ --include-before-body='before_body.html'  $^

watch:
	watchman watch .
	watchman -- trigger . build-html Makefile '*.md' before_body.html -- make html
	watchman -- trigger . build-pdf Makefile '*.tex' -- make pdf

unwatch:
	watchman watch-del .
