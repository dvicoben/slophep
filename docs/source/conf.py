import sys
import os
sys.path.insert(0, os.path.abspath('../../src/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SLOP'
copyright = '2025, dvicoben'
author = 'dvicoben'
release = 'v1.3.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',            # To generate autodocs
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',           # autodoc with maths
    'sphinx.ext.napoleon',           # For auto-doc configuration
    'sphinx.ext.viewcode'
]

templates_path = ['_templates']
exclude_patterns = []
autosummary_generate=True
# If true, the current module name will be prepended to all description unit titles (such as .. function::).
# add_module_names = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']
