.PHONY: all

TEX=$(shell find tex -name '*.tex')
TXT=$(shell find tex -name '*.txt')
PDF=$(subst .tex,.pdf,$(shell find tex -name '*.tex' | grep -v preamble))

.PHONY: all
all: $(PDF)

$(PDF): $(TEX) $(TXT)
	$(MAKE) -C tex $(notdir $@)
