PY = $(shell find ../python -name '*.py')
TOP = $(addprefix ../,flake.nix flake.lock Makefile env.sh)
LIB = lib/preamble.tex lib/template.tex lib/references.bib

%.tar.gz: %.tex %.yaml Makefile $(LIB) $(PY) $(TOP)
	mksource.sh $<

%.pdf: %.tex %.yaml %.tar.gz Makefile $(LIB) $(PY) $(TOP)
	gpg --armor --yes --output $*.tar.gz.sig --detach-sign $*.tar.gz
	sed -e 's@!BODYFILE!@$*.tex@g' \
	    -e 's@!TITLE!@$(shell yq -r '.title' $*.yaml)@g' \
	    -e 's@!ABSTRACT!@$(shell yq '.abstract' $*.yaml | sed 's@^.@@;s@.$$@@')@g' \
	    -e 's@!GITHASH!@$(shell git rev-parse --short=7 HEAD)@g' \
	    -e 's@!GPGSIG!@$*.tar.gz.sig@g' \
	    -e 's*!EMAIL!*$(shell gpg --list-keys --with-colons | grep "^uid" | head -n 1 | awk -F "<|>" "{print \$$2}")*g' \
	    -e 's@!SOURCELINK!@$(shell texlink.py $<)@g' \
		lib/template.tex > _note.tex
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error -outdir=_build  _note.tex
	cp _build/_note.pdf $*.pdf

# Does not work properly
%.signed.pdf: %.pdf
	gpg --clearsign --output=$@ $<

