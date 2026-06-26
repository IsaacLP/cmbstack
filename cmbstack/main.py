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
    """

    def __init__(self, sky_map, nside):
        self.map = sky_map
        self.nside = nside
        self.normalized = None
        self.positions = None
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
        sky_map = np.asarray(sky_map)
        nside = hp.npix2nside(sky_map.size)
        return cls(sky_map, nside)
 
    @classmethod
    def from_fits(cls, path, field=0):
        """Build a pipeline from any HEALPix FITS file on disk."""
        sky_map = maps.load_map(path, field=field)
        return cls.from_map(sky_map)
    

    def run(self, size_deg=10.0, reso_arcmin=3.0, profile=True, threshold=3.0, n_peaks=None):
        """Run the full stacking loop.

        Parameters
        ----------
        size_deg, reso_arcmin : float
            Patch geometry.
        profile : bool
            Whether to also compute the radial profile.
        threshold : float
            Peak-finding threshold in units of the map std.
        n_peaks : int or None
            Maximum number of peaks to use; None means use all.

        Returns
        -------
        result : cmbstack.stack.StackResult
        """

        self.normalized = maps.normalize_map(self.map)

        self.positions = stacking.find_peaks(self.normalized, self.nside, threshold=threshold, n_peaks=n_peaks)

        self.patches = stacking.extract_patches(self.normalized,self.positions,size_deg=size_deg,reso_arcmin=reso_arcmin)

        self.stacked = stacking.stack_patches(self.patches)

        if profile:
            self.radius, self.profile = stacking.radial_profile(self.stacked,reso_arcmin=reso_arcmin)

        return self