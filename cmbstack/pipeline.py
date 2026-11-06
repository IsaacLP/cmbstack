"""High-level stacking pipeline for HEALPix CMB maps.

:class:`StackingPipeline` runs the full analysis in one place:
simulate (or load) a map, detect peaks, extract gnomonic patches, stack them,
and compute a radial profile. The pipeline stores every intermediate product
as an attribute so individual steps can be inspected after the run.

Typical use
-----------
>>> pipeline = StackingPipeline.from_cl("path/to/cl.txt", nside=1024, seed=42)
>>> pipeline.run()
>>> plt.imshow(pipeline.stacked)
"""

from . import maps, stacking
import numpy as np
import healpy as hp


class StackingPipeline:
    """End-to-end stacking pipeline.

    Construct from a power spectrum (:meth:`from_cl`) or from an existing map
    (:meth:`from_map`), then call :meth:`run`.

    Parameters
    ----------
    sky_map : numpy.ndarray
        The HEALPix map to stack on.
    nside : int
        Resolution parameter of the map.

    Attributes
    ----------
    normalized : numpy.ndarray or None
        Set after run(); the normalized map.
    positions : numpy.ndarray or None
        Set after run(); the peak or catalogue positions in (lon, lat) degrees.
    psi : numpy.ndarray or None
        Set after run(); the per-patch rotation angles in degrees, if any.
    patches : numpy.ndarray or None
        Set after run(); the extracted patches.
    stacked : numpy.ndarray or None
        Set after run(); the stacked image.
    radius : numpy.ndarray or None
        Set after run(); the radial distances of the profile in arcmin.
    profile : numpy.ndarray or None
        Set after run(); the radial profile of the stacked image.
    """

    def __init__(self, sky_map, nside):
        self.map = sky_map
        self.nside = nside
        self.normalized = None
        self.positions = None
        self.psi = None
        self.patches = None
        self.stacked = None
        self.radius = None
        self.profile = None


    @classmethod
    def from_cl(cls, cl_path, nside=128, seed=None):
        """Build a pipeline by simulating a map from a power-spectrum file."""
        cl = maps.load_cl(cl_path)

        sky_map = maps.simulate_map(cl, nside, seed)
        return cls(sky_map, nside)
    
    @classmethod
    def from_map(cls, sky_map):
        """Build a pipeline from a HEALPix map array already in memory. nside is inferred from the map length, so the caller doesn't have to supply it.
        """
        # Masked pixels become UNSEEN; np.asarray alone would drop the mask.
        sky_map = np.ma.filled(sky_map, hp.UNSEEN)
        nside = hp.npix2nside(sky_map.size)
        return cls(sky_map, nside)
 
    @classmethod
    def from_fits(cls, path, field=0):
        """Build a pipeline from any HEALPix FITS file on disk."""
        sky_map = maps.load_map(path, field=field)
        return cls.from_map(sky_map)
    

    def run(self, size_deg=10.0, reso_arcmin=3.0, profile=True, threshold=3.0, n_peaks=None,
            positions=None, psi=None, min_count=1, mask_distance_deg=None, verbose=False,
            size_deg_y=None, keep_patches=True):
        """Run the full stacking loop.

        Parameters
        ----------
        size_deg, reso_arcmin : float
            Patch geometry: side length in x (degrees) and pixel size (arcmin).
        size_deg_y : float or None
            Patch side length in y in degrees. None (default) gives a square
            patch of side ``size_deg``.
        profile : bool
            Whether to also compute the radial profile.
        threshold : float
            Peak-finding threshold in units of the map std.
        n_peaks : int or None
            Maximum number of peaks to use; None means use all.
        positions : array_like, shape (N, 2), optional
            External catalogue positions as ``(lon, lat)`` in degrees. 
            If given, peak finding is skipped and ``threshold``/``n_peaks`` 
            are ignored.
        psi : array_like, shape (N,), optional
            Per-patch rotation in degrees for oriented stacking.
        min_count : int or float
            Minimum number of contributing patches for a stacked pixel to be
            kept (see :func:`stacking.stack_patches`). The per-pixel counts
            are stored in ``self.counts``.
        mask_distance_deg : float or None
            Discard detected peaks with a masked pixel within this distance in
            degrees (see :func:`stacking.find_peaks`). Peaks touching the mask
            are always discarded. Ignored when ``positions`` is given. The
            number of discarded peaks is stored in ``self.n_rejected``.
        verbose : bool
            If True, print how many peaks were kept and rejected by the mask cut.
        keep_patches : bool
            If True (default), keep every patch in ``self.patches``. If False,
            stream the patches straight into the stack and set ``self.patches``
            to None: memory then stays at about one patch whatever the number
            of positions, which is needed for large catalogues.

        Returns
        -------
        self : StackingPipeline
            The pipeline object itself, with all intermediate products stored as
            attributes.
        """

        self.normalized = maps.normalize_map(self.map)

        if positions is None:
            self.positions, self.n_rejected = stacking.find_peaks(
                self.normalized, self.nside, threshold=threshold, n_peaks=n_peaks,
                mask_distance_deg=mask_distance_deg, return_n_rejected=True)
            if verbose:
                print(f"find_peaks: kept {len(self.positions)} peaks, "
                      f"rejected {self.n_rejected} near the mask")
        else:
            self.positions = np.atleast_2d(positions)
            self.n_rejected = 0
        if len(self.positions) == 0:
            raise ValueError("No positions to stack on: no peaks found above the threshold.")
        self.psi = psi

        patches = stacking.iter_patches(self.normalized,self.positions,psi=self.psi,size_deg=size_deg,
                                        size_deg_y=size_deg_y,reso_arcmin=reso_arcmin)
        if keep_patches:
            patches = list(patches)
        self.patches = patches if keep_patches else None

        self.stacked, self.counts = stacking.stack_patches(patches, min_count=min_count, return_counts=True)

        if profile:
            self.radius, self.profile = stacking.radial_profile(self.stacked,reso_arcmin=reso_arcmin)

        return self

StackingPipeline.__module__ = "cmbstack"