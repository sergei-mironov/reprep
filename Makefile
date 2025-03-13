.PHONY: all

TEX = $(shell find tex -name '[^_]*.tex')
YAML = $(shell find tex -name '*.yaml')
PY = $(shell find $(shell pwd)/python -name '*.py')
PDF = $(subst .yaml,.pdf,$(YAML))
DEP = $(subst .pdf,.pdf.dep,$(PDF))
TOP = flake.nix flake.lock Makefile env.sh

.PHONY: all
all: $(PDF)

$(PDF): $(TEX) $(TXT) $(YAML) $(PY) $(TOP)
  $(MAKE) -C tex $(notdir $@)

$(DEP):
	-$(MAKE) -p -C tex $(notdir $@)

.PHONY:help
help:
	echo $(PY)
	echo $(DEP)
	$(MAKE) -C tex help
