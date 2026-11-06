cmbstack
====================================

cmbstack is a Python package for stacking patches of `HEALPix <https://healpy.readthedocs.io/en/latest/index.html>`_ maps, such as the Cosmic Microwave Background (CMB) temperature, CMB lensing convergence or galaxy density. It accepts as input a theoretical power spectrum, a HEALPix FITS file, or a map array already in memory, which can be masked. It stacks either on the local maxima of the map or on the positions of an external catalogue, with an optional rotation for each patch. It extracts gnomonic (flat-sky) patches around each position and averages them. This stacking procedure enhances the coherent profile while suppressing uncorrelated noise.

The tutorials below are worked examples. The notebooks and the files required to run them can be downloaded `here <https://github.com/IsaacLP/cmbstack/tree/main/docs/tutorials>`_:

- **Simulated map:** every step of the pipeline on a simulated CMB temperature map.
- **Planck peak stacking:** peak stacking on the masked Planck SMICA map.
- **Oriented stacking on filaments:** stacking an SDSS galaxy overdensity map on a cosmic filament catalogue, with each patch rotated along its filament.


.. toctree::
   :maxdepth: 1
   :caption: User guide

   installation
   Tutorial: simulated map <tutorials/tutorial>
   Tutorial: Planck peak stacking <tutorials/planck_peak_stacking>
   Tutorial: oriented stacking on filaments <tutorials/filament_oriented_stacking>
   api
