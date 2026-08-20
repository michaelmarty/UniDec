# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import tomllib
from datetime import date
from pathlib import Path

# Make the public repository importable on Windows and in the Linux Pages runner.
repository_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repository_root))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with (repository_root / 'pyproject.toml').open('rb') as pyproject_file:
    project_metadata = tomllib.load(pyproject_file)['project']

project = project_metadata['name']
copyright = f'{date.today().year}, University of Arizona'
author = 'Michael Marty'
release = project_metadata['version']
version = release

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
# Run unidec_doc/generate.bat on Windows to regenerate the API pages and build HTML.

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'myst_parser']

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
    'unidec.IsoDec.IsoGen.isogenc',
    'unidec.IsoDec.isogenwrapper',
    'unidec.modules.unidecwrapper',
    'wx',
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
