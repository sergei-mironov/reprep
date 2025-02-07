# %.pdf: %.tex preamble.tex
# 	latexmk -shell-escape -pdf -interaction=nonstopmode \
# 	        -latex=pdflatex --halt-on-error -outdir=_build $*.tex
# 	cp _build/$*.pdf .

%.pdf: %.tex %.yaml preamble.tex note-template.tex Makefile
	sed -e 's/!BODY!/$*.tex/g' \
	    -e 's/!TITLE!/$(shell yq -r '.title' $*.yaml)/g' \
	    -e 's/!ABSTRACT!/$(shell yq '.abstract' $*.yaml | sed 's/^.//;s/.$$//')/g' note-template.tex > _note.tex
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error -outdir=_build  _note.tex
	cp _build/_note.pdf $*.pdf

