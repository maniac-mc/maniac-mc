# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MANIAC-MC'
copyright = '2025, Simon Gravelle'
author = 'Simon Gravelle'
release = 'v0.1.2-alpha'

import os
breathe_projects = {
    "ManiaC": os.path.abspath(os.path.join(os.path.dirname(__file__), "xml"))
}
print("Breathe XML path:", breathe_projects["ManiaC"])
breathe_default_project = "ManiaC"

# -- General configuration ---------------------------------------------------

extensions = [
    'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'
html_static_path = ['_static']
