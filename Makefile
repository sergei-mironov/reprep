.PHONY: all

all: tex/main.pdf

tex/%.pdf: tex/preamble.tex tex/%.tex
	cd tex && \
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error $*.tex

_ai/%.pdf: tex/preamble.tex _ai/%.tex
	cd _ai && \
	latexmk -shell-escape -pdf -interaction=nonstopmode \
	        -latex=pdflatex --halt-on-error $*.tex
