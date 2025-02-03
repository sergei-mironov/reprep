.PHONY: all

TEX=$(shell find tex -name '*.tex')
PDF=$(subst .tex,.pdf,$(shell find tex -name '*.tex' | grep -v preamble))

.PHONY: all
all: $(PDF)

$(PDF): $(TEX)
	$(MAKE) -C tex $(notdir $@)
