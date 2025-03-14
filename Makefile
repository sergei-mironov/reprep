.PHONY: all

TEX = $(shell find tex -name '[^_]*.tex')
YAML = $(shell find tex -name '*.yaml')
PDF = $(subst .yaml,.pdf,$(YAML))

.PHONY: all $(PDF)
all: $(PDF)

$(PDF):
	$(MAKE) -C tex $(notdir $@)

$(DEP):
	-$(MAKE) -p -C tex $(notdir $@)

.PHONY:help
help:
	$(MAKE) -C tex help
