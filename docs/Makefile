# Docs Makefile - builds Sphinx docs only

# Sphinx configuration
SPHINXBUILD   ?= sphinx-build
SPHINXOPTS    ?=
SOURCEDIR     = source
BUILDDIR      = build

.PHONY: help html clean all

# Default target: build HTML docs
all: html

# Show help
help:
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  html      - build HTML docs"
	@echo "  clean     - remove generated files"

# Build HTML docs
html:
	@echo "Building Sphinx HTML..."
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)
	@echo "HTML docs built in $(BUILDDIR)/html"

# Clean Sphinx build
clean:
	@echo "Cleaning Sphinx build..."
	@rm -rf $(BUILDDIR)
	@echo "Clean complete"
