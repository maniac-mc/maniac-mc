# Docs Makefile - builds Doxygen + Sphinx docs

# Sphinx configuration
SPHINXBUILD   ?= sphinx-build
SPHINXOPTS    ?=
SOURCEDIR     = source
BUILDDIR      = build

# Doxygen configuration
DOXYFILE      ?= doxyfile
DOXYGEN_XML   = $(SOURCEDIR)/xml

.PHONY: help doxygen html clean all

# Default target: build HTML docs
all: html

# Show help
help:
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  html      - build HTML docs (runs Doxygen first)"
	@echo "  doxygen   - generate Doxygen XML for Sphinx"
	@echo "  clean     - remove generated files"

# Generate Doxygen XML
doxygen:
	@echo "Generating Doxygen XML..."
	@doxygen $(DOXYFILE)

# Build HTML docs (runs Doxygen first)
html: doxygen
	@echo "Building Sphinx HTML..."
	@$(SPHINXBUILD) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS)
	@echo "HTML docs built in $(BUILDDIR)/html"

# Clean Sphinx build and Doxygen XML
clean:
	@echo "Cleaning Sphinx build and Doxygen XML..."
	@rm -rf $(BUILDDIR) $(DOXYGEN_XML)
	@echo "Clean complete"
