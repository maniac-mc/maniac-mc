#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# MANIAC-MC Documentation Build Script
# ------------------------------------------------------------------------------
# 1. Reads the version number from version.txt (repository root)
# 2. Updates the 'release' variable in Sphinx conf.py
# 3. Builds the documentation using the Makefile (Doxygen + Sphinx)
# ------------------------------------------------------------------------------

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DOCS_DIR="$ROOT_DIR/docs"
CONF_FILE="$DOCS_DIR/source/conf.py"
VERSION_FILE="$ROOT_DIR/version.txt"

# --- Step 1: Read version number ---------------------------------------------
if [[ ! -f "$VERSION_FILE" ]]; then
    echo "Error: version.txt not found at $VERSION_FILE"
    exit 1
fi

VERSION=$(tr -d '[:space:]' < "$VERSION_FILE")
echo "Using version: $VERSION"

# --- Step 2: Update version in conf.py ---------------------------------------
if [[ ! -f "$CONF_FILE" ]]; then
    echo "Error: conf.py not found at $CONF_FILE"
    exit 1
fi

echo "Updating conf.py with version: $VERSION"
sed -i "s/^release = .*/release = '${VERSION}'/" "$CONF_FILE"

# --- Step 3: Build documentation ---------------------------------------------
cd "$DOCS_DIR"
echo "Running 'make html'..."
make html

echo "Documentation build complete."
