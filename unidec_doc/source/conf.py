# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import tomllib
from pathlib import Path

# Make the repository importable on Windows and in the Linux Pages runner.
repository_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repository_root))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with (repository_root / 'pyproject.toml').open('rb') as pyproject_file:
    project_metadata = tomllib.load(pyproject_file)['project']

project = project_metadata['name']
copyright = '2026, University of Arizona'
author = 'Michael Marty'
release = project_metadata['version']
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# cd to unidec_doc
# sphinx-apidoc -f -o source ..\unidec
# .\make.bat html

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
