PDF_FILES = ssmodels_in_stan.pdf
HTML_FILES = ssmodels.html ssfunctions.md
BIB = ../bib/jrnold.bib

build: pdf html

pdf: $(PDF_FILE)

html: $(HTML_FILES)

%.pdf: %.tex
	latexmk -pdf -xelatex $^

%.html: %.md
	pandoc -s --mathjax -f markdown+pipe_tables -o $@ --include-before-body='before_body.html'  $^

watch:
	watchman watch .
	watchman -- trigger . build-html Makefile '*.md' before_body.html -- make html
	watchman -- trigger . build-pdf Makefile '*.tex' -- make pdf

unwatch:
	watchman watch-del .
