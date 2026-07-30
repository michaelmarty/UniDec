# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'UniDec'
copyright = '2025, University of Arizona'
author = 'Michael Marty'
release = '8.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# cd to unidec_doc
# sphinx-apidoc -f -o source ..\unidec
# .\make.bat html

import sys
from pathlib import Path

# Make the repository importable on Windows and in the Linux Pages runner.
repository_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repository_root))

#for x in os.walk('..\\..\\unidec'):
  #sys.path.insert(0, x[0])

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.coverage', 'myst_parser']

# Avoid importing optional GUI, training, vendor, and platform-specific
# dependencies on the documentation runner.
autodoc_mock_imports = [
    'Scripts',
    'System',
    'ThermoFisher',
    'clr',
    'isogen_tools',
    'isogenatom',
    'isogenatom_trainingdata',
    'isogenpep_trainingdata',
    'massql',
    'matchms',
    'molmass',
    'mpld3',
    'sqlalchemy',
    'torch',
    'torchvision',
    'wx',
]

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
