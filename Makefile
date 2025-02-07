.PHONY: all

TEX=$(shell find tex -name '[^_]*.tex')
YAML=$(shell find tex -name '*.yaml')
PDF=$(subst .yaml,.pdf,$(YAML))

.PHONY: all
all: $(PDF)

$(PDF): $(TEX) $(TXT) $(YAML)
	$(MAKE) -C tex $(notdir $@)
