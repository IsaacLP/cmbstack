.. cmbstack documentation master file, created by
   sphinx-quickstart on Fri Jun 26 00:23:37 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

"cmbstack" documentation:
===========================

Welcome to cmbstack's documentation! 
====================================

This **cmbstack** is a package built with Python programming language, for stacking the Cosmic microwave Background (CMB) maps around data around galaxy clusters data - so that we can then study the thermal Sunyaev-Zel'dovich (i.e. the tSZ) effect using this stacking result.

Key features of this cmbstack package:
- This package generates the simulated CMB temperature maps using the user's given input power spectra.
- While working, this package stacks the generated simulated CMB temperature maps around the galaxy cluster centers.
- This package visualizes the results obtained from stacking with customisable plots.


Installation Procedure:
-----------------------

Option-[1]: 
Install the latest version from PyPI:

.. code-block:: bash

   pip install cmbstack

Or 

Option-[2]:
Install the development version from GitHub:

.. code-block:: bash

   git clone https://github.com/IsaacLP/cmbstack.git
   cd cmbstack
   pip install -e .


Quick Start Method:
-------------------

.. code-block:: python

   import cmbstack as cs
   # Your code here

For more detailed usage, please see the API Reference section below.



Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api 

