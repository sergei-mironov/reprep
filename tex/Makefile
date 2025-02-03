%.pdf: %.tex preamble.tex
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error -outdir=_build $*.tex
	cp _build/$*.pdf .

