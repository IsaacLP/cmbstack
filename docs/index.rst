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


Changelog
---------

v0.2.0
~~~~~~

**Breaking changes**

- Sky positions are now ``(lon, lat)`` in degrees everywhere (``find_peaks``, ``extract_patches``, ``StackingPipeline.run``), instead of ``(theta, phi)`` in radians.

**New features**

- ``StackingPipeline.run`` accepts external ``positions`` (previously only possible by calling ``extract_patches`` directly) and per-patch rotation angles ``psi``, for oriented stacking.
- Rectangular patches: ``size_deg_y`` in ``extract_patches`` and ``run()``, with ``size_deg`` as the side along x.
- ``size_mpc_to_deg`` converts a physical size in Mpc to an angle at redshift z, for an astropy cosmology.
- Masked maps:

  - ``from_map`` keeps the mask of an ``hp.ma`` array.
  - ``normalize_map`` computes the mean and std on unmasked pixels only, including for float32 maps.
  - Masked pixels become NaN in the patches.

- Per-pixel hit counts: ``stack_patches`` divides each pixel by the number of patches with a valid value there, takes ``min_count``, and can return the count map (``return_counts=True``, ``pipeline.counts``).
- ``find_peaks`` drops peaks next to the mask, and ``mask_distance_deg`` widens the cut. ``return_n_rejected=True``, ``run(verbose=True)`` and ``pipeline.n_rejected`` report how many were dropped.
- Streaming: ``iter_patches`` yields one patch at a time, ``stack_patches`` accepts any iterable, and ``run(keep_patches=False)`` stacks any number of positions in constant memory.

**Other changes**

- ``extract_patches`` uses ``GnomonicProj`` directly, which is faster because it no longer creates a matplotlib figure per patch.
- ``run()`` raises a clear ``ValueError`` when no peaks are found.
- New tutorials: peak stacking on the Planck SMICA map, and oriented stacking of an SDSS galaxy overdensity map on a cosmic filament catalogue.

v0.1.1
~~~~~~

- ``load_cl`` reads the chosen column with pandas.
- Runtime test placeholder moved to ``tests/test_runtime.py``.
