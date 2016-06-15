PDF_FILE = ssmodels_in_stan.pdf
HTML_FILE = ssmodels.html

build: pdf html

pdf: $(PDF_FILE)

html: $(HTML_FILE)

$(PDF_FILE): $(PDF_FILE:%.pdf=%.tex)
	latexmk -pdf -xelatex $^

$(HTML_FILE): $(HTML_FILE:%.html=%.md)
	pandoc -S -o $@ $^

watch:
	watchman watch .
	watchman -- trigger . build '*.tex' '*.html' -- make build

unwatch:
	watchman watch-del .
