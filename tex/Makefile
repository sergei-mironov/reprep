PY ?= $(shell find ../python -name '*.py')

%.pdf: %.tex %.yaml lib/preamble.tex lib/template.tex lib/references.bib Makefile $(PY)
	sed -e 's@!BODYFILE!@$*.tex@g' \
	    -e 's@!TITLE!@$(shell yq -r '.title' $*.yaml)@g' \
	    -e 's@!ABSTRACT!@$(shell yq '.abstract' $*.yaml | sed 's@^.@@;s@.$$@@')@g' \
	    -e 's@!GITHASH!@$(shell git rev-parse --short=7 HEAD)@g' \
	    -e 's@!SOURCELINK!@$(shell texlink.py $<)@g' \
			lib/template.tex > _note.tex
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error -outdir=_build  _note.tex
	cp _build/_note.pdf $*.pdf
