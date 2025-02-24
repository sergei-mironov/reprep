.PHONY: all

TEX=$(shell find tex -name '[^_]*.tex')
YAML=$(shell find tex -name '*.yaml')
PY=$(shell find $(shell pwd)/python -name '*.py')
PDF=$(subst .yaml,.pdf,$(YAML))

export PY

.PHONY: all
all: $(PDF)

$(PDF): $(TEX) $(TXT) $(YAML) $(PY)
	$(MAKE) -C tex $(notdir $@)

.PHONY:help
help:
	echo $(PY)
	$(MAKE) -C tex help
