# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'UniDec'
copyright = '2025, University of Arizona'
author = 'Michael Marty'
release = '8.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# cd to unidec_doc
# sphinx-apidoc -f -o source ..\unidec
# .\make.bat html

import sys, os

sys.path.append("C:\\Python\\UniDec3\\unidec")

topdir = os.path.abspath('..\\..\\unidec')
print(topdir)
sys.path.insert(0, topdir)

topdir = os.path.abspath('..\\..')
print(topdir)
sys.path.insert(0, topdir)

#for x in os.walk('..\\..\\unidec'):
  #sys.path.insert(0, x[0])

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.coverage']  # , 'myst_parser']

templates_path = ['_templates']
# exclude_patterns = ["**/thermo_reader", "**/waters_importer"]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_static_path = ['_static']
