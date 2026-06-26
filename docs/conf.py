# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import sys

# ----- Add your package path here below ----------
# This will tell Sphinx where to find 'cmbstack' python code.
# '..' means 'go up one folder' == this will take Sphinx from the current docs folder to the main project root.
sys.path.insert(0, os.path.abspath('..'))



# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cmbstack'
copyright = '2026, Isaac Alexis López Paredes, Anushka Sanjay Tilekar'
author = 'Isaac Alexis López Paredes, Anushka Sanjay Tilekar'
release = '0.0.1'
root_doc = 'index'  # root_doc tells Sphinx which .rst file is your main homepage. 
                    # It's set to 'index' by default (which matches your index.rst file), 
                    # so you're fine either way. If you don't add it, 
                    # Sphinx will just use the default.




# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',    # This extention Automatically reads docstrings from our code
    'sphinx.ext.napoleon',   # This extensoin allows us to use Google style docstrings
    'sphinx.ext.mathjax',    # This extension renders the mathematical equations - because we have used a lot of mathematical equations in this repository, haha!
    'sphinx.ext.viewcode',   # This extension will add links to our source code on Github
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
# Trying a different looking ReadTheDocs theme instead of the default 'alabaster' taught in lecture
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
