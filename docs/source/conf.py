import sys
import os
sys.path.insert(0, os.path.abspath('../../src/'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project   = 'SLOP'
copyright = '2026, dvicoben'
author    = 'dvicoben'
release   = 'v1.5.1'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "myst_parser",
    'sphinx.ext.autodoc',            # To generate autodocs
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',           # autodoc with maths
    'sphinx.ext.napoleon',           # For auto-doc configuration
    'sphinx.ext.viewcode'
]

myst_enable_extensions = [
    "dollarmath", 
    "amsmath"
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md' : 'markdown',
}

templates_path       = ['_templates']
exclude_patterns     = []
autosummary_generate = True
# If true, the current module name will be prepended to all description unit titles (such as .. function::).
add_module_names = False


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
# html_static_path = ['_static']
html_context = {
    # Show the "Edit on Git" link instead of "View Source"
    "display_github" : True,
    "github_host"    : "github.com",
    "github_user"    : "dvicoben",
    "github_repo"    : "slophep",
    "github_version" : "master/docs/source/",
}
rst_prolog = """
:github_url: https://github.com/dvicoben/slophep
"""
