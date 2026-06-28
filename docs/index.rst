cmbstack 
====================================

cmbstack is a Python package for stacking patches of the Cosmic Microwave Background (CMB) temperature sky. It accepts as input a theoretical power spectrum, a `HEALPix <https://healpy.readthedocs.io/en/latest/index.html>`_ FITS file, or a map array already in memory. From there it detects local maxima, extracts gnomonic (flat-sky) patches around each peak, and averages them. This stacking procedure enhances the coherent peak profile while suppressing uncorrelated noise.

For a full worked example, see the Tutorial page below or download the jupyter notebook version and the files required to run it `here <https://github.com/IsaacLP/cmbstack/tree/main/docs/tutorials>`_.


.. toctree::
   :maxdepth: 1
   :caption: User guide

   installation
   Tutorial <tutorials/tutorial>
   api