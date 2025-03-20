R = ..
PY = $(shell find ../python -name '*.py')
SH = $(shell find ../sh -name '*')
TOP = $(addprefix ../,flake.nix flake.lock Makefile env.sh)
LIB = $(shell find ./lib -name '*tex')

%.tar.gz: %.tex %.yaml Makefile $(LIB) $(PY) $(TOP) $(SH)
	cd $R && tar -czf \
		$(shell realpath --relative-to=$R $@) \
	  $(foreach prereq,$^,$(shell realpath --relative-to=$R $(prereq)))

%.sig: %.tar.gz
	gpg --armor --yes --output $@ --detach-sign $<

%.pdf: %.sig
	sed -e 's@!BODYFILE!@$*.tex@g' \
	    -e 's@!DESTFILE!@$@@g' \
	    -e 's@!TITLE!@$(shell yq -r '.title' $*.yaml)@g' \
	    -e 's@!ABSTRACT!@$(shell yq '.abstract' $*.yaml | sed 's@^.@@;s@.$$@@')@g' \
	    -e 's@!GITHASH!@$(shell git rev-parse --short=7 HEAD)@g' \
	    -e 's@!GPGSIGFILE!@$*.sig@g' \
	    -e 's*!GPGNAME!*$(shell gpgdata.sh | head -n1 | tail -n1)*g' \
	    -e 's*!GPGEMAIL!*$(shell gpgdata.sh | head -n2 | tail -n1)*g' \
	    -e 's@!SOURCELINK!@$(shell texlink.py $<)@g' \
	    -e 's@!SOURCEFILE!@$*.tar.gz@g' \
		lib/template.tex > _note.tex
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error -outdir=_build  _note.tex
	cp _build/_note.pdf $*.pdf

# Does not work properly
%.signed.pdf: %.pdf
	gpg --clearsign --output=$@ $<

