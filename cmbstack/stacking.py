"""Stacking of patches around positions on a HEALPix map.

The pipeline is field-agnostic: it operates on any scalar HEALPix map (CMB
temperature, lensing convergence, a y-map, galaxy density, ...). Positions to
stack on can either be auto-detected peaks (local maxima) or supplied as an
external catalogue, so the same code serves peak stacking and
stacking-on-catalogue (clusters, voids, filaments, ...).

Typical use
-----------
>>> peaks = find_peaks(m, nside, threshold=3.0)      # positions in (theta, phi)
>>> patches = extract_patches(m, peaks)              # fixed-grid 2D cutouts
>>> stacked = stack_patches(patches)                 # mean 2D image
>>> r, profile = radial_profile(stacked, reso_arcmin=3.0)
"""

import numpy as np
import healpy as hp


# ---------------------------------------------------------------------------
# Position finding
# ---------------------------------------------------------------------------
def find_peaks(sky_map, nside, threshold=None, n_peaks=None, mask_distance_deg=None,
               return_n_rejected=False):
    """Find local maxima of a HEALPix map and return their sky positions.

    A pixel is a local maximum if its value is strictly greater than all of its
    immediate HEALPix neighbours. Peaks can be filtered by a significance
    threshold and/or capped at the ``n_peaks`` highest.

    Peaks next to masked (UNSEEN) pixels are discarded: there a pixel only has
    to beat its unmasked neighbours, so the true maximum may lie under the mask.

    Parameters
    ----------
    sky_map : numpy.ndarray
        Input HEALPix map (RING ordering). For a normalised map, values are in
        units of sigma, so ``threshold`` is a significance nu.
    nside : int
        HEALPix resolution parameter of ``sky_map``.
    threshold : float, optional
        If given, keep only peaks with value greater than this (e.g. 3.0 for
        3-sigma peaks on a normalised map). Default None (no threshold).
    n_peaks : int, optional
        If given, keep only the ``n_peaks`` highest peaks (applied after the
        threshold and the mask cut). Default None (keep all).
    mask_distance_deg : float, optional
        If given, also discard peaks with any masked pixel within this angular
        distance in degrees. Default None (only the immediate neighbours are
        checked).
    return_n_rejected : bool, optional
        If True, also return the number of peaks discarded for being next to
        (or within ``mask_distance_deg`` of) the mask. Default False.

    Returns
    -------
    positions : numpy.ndarray, shape (N, 2)
        Sky positions of the selected peaks as ``(lon, lat)`` in degrees, the
        same format accepted by :func:`extract_patches`.
    n_rejected : int
        Number of peaks discarded by the mask cut (before ``threshold`` and
        ``n_peaks``). Only returned if ``return_n_rejected`` is True.
    """
    # hp.hotspots returns (max_map, minima_pix, maxima_pix); we want the maxima.
    _, _, maxima_pix = hp.hotspots(sky_map)
    maxima_pix = np.asarray(maxima_pix)

    # Drop peaks touching the mask (or within mask_distance_deg of it).
    bad = hp.mask_bad(sky_map)
    n_rejected = 0
    if bad.any() and maxima_pix.size:
        neighbours = hp.get_all_neighbours(nside, maxima_pix)
        near_mask = np.any((neighbours >= 0) & bad[neighbours], axis=0)
        if mask_distance_deg is not None:
            radius = np.radians(mask_distance_deg)
            vecs = hp.pix2vec(nside, maxima_pix)
            near_mask |= [
                bad[hp.query_disc(nside, v, radius, inclusive=True)].any()
                for v in np.transpose(vecs)
            ]
        n_rejected = int(np.count_nonzero(near_mask))
        maxima_pix = maxima_pix[~near_mask]

    values = sky_map[maxima_pix]

    if threshold is not None: # Add check for minimmun threshold/n_peaks and ratio between threshold and n_peaks
        keep = values > threshold
        maxima_pix = maxima_pix[keep]
        values = values[keep]

    # Sort highest-first so n_peaks keeps the most significant.
    order = np.argsort(values)[::-1]
    maxima_pix = maxima_pix[order]

    if n_peaks is not None:
        maxima_pix = maxima_pix[:n_peaks]

    lon, lat = hp.pix2ang(nside, maxima_pix,lonlat=True)
    positions = np.column_stack([lon, lat])
    if return_n_rejected:
        return positions, n_rejected
    return positions


# ---------------------------------------------------------------------------
# Patch extraction
# ---------------------------------------------------------------------------
def size_mpc_to_deg(size_mpc, z, cosmo):
    """Convert a physical size in Mpc to an angular size in degrees.

    Uses the angular diameter distance at redshift ``z`` from the given
    cosmology.

    Parameters
    ----------
    size_mpc : float
        Physical size in Mpc.
    z : float
        Redshift.
    cosmo : astropy.cosmology.Cosmology
        Cosmology object to compute the angular diameter distance.

    Returns
    -------
    size_deg : float
        Angular size in degrees.
    """
    d_a = cosmo.angular_diameter_distance(z).value  # Mpc
    size_rad = size_mpc / d_a  # radians
    size_deg = np.degrees(size_rad)
    return size_deg


def extract_patches(sky_map, positions, psi=None, size_deg_x=10.0, size_deg_y=None, reso_arcmin=3.0):
    """Extract fixed-grid gnomonic patches centred on each position.

    Each patch is a square or a rectangular 2D array produced by a gnomonic (tangent-plane)
    projection centred on the position, so every patch shares the same grid and
    the centre pixel always corresponds to the position itself. 

    Parameters
    ----------
    sky_map : numpy.ndarray
        Input HEALPix map (RING ordering). Any scalar field.
    positions : array_like, shape (N, 2)
        Sky positions as ``(lon, lat)`` in degrees (e.g. the output of
        :func:`find_peaks`, or an external catalogue converted to this format).
    psi : array_like, shape (N,), optional
        Rotation angle in degrees for each patch, passed to the gnomonic projector's
        ``rot``. If None, defaults to 0.0.
    size_deg_x : float, optional
        Full side length of the square patch in degrees. Default 10.0.
    size_deg_y : float, optional
        Full side length of the rectangular patch in degrees. If None, defaults
        to ``size_deg_x`` (square patch). Default None.
    reso_arcmin : float, optional
        Pixel size of the projected patch in arcminutes. Default 3.0.

    Returns
    -------
    patches : list of numpy.ndarray
        One square or rectangular 2D array per position, all of identical shape
        if 
        ``(xsize, xsize)`` with ``xsize = size_deg * 60 / reso_arcmin``.
    """
    positions = np.atleast_2d(positions)
    if size_deg_y is None:
        size_deg_y = size_deg_x
    xsize = int(round(size_deg_x * 60.0 / reso_arcmin))
    ysize = int(round(size_deg_y * 60.0 / reso_arcmin))

    # Masked arrays become plain arrays with UNSEEN at the masked pixels.
    sky_map = hp.pixelfunc.ma_to_array(sky_map)
    nside = hp.get_nside(sky_map)

    def vec2pix(x, y, z):
        return hp.vec2pix(nside, x, y, z)

    patches = []
    if psi is None:
        psi = np.zeros(len(positions))
    for lon, lat, angle in zip(positions[:, 0], positions[:, 1], psi):
        # Use the projector directly: same result as gnomview, but without
        # building a matplotlib figure per patch. rot expects (lon, lat) in degrees.
        proj = hp.projector.GnomonicProj(
            rot=(lon, lat, angle), xsize=xsize, ysize=ysize, reso=reso_arcmin
        )
        patch = proj.projmap(sky_map, vec2pix).astype(float)
        # Masked (UNSEEN) pixels become NaN so stack_patches ignores them.
        patch[hp.mask_bad(patch)] = np.nan
        patches.append(patch)

    return patches


# ---------------------------------------------------------------------------
# Stacking
# ---------------------------------------------------------------------------
def stack_patches(patches, min_count=1, return_counts=False):
    """Average patches pixel-by-pixel into a single stacked image.

    Because every patch shares the same fixed grid (see :func:`extract_patches`),
    the mean is a genuine stacked image: incoherent noise averages towards zero
    while the coherent central profile survives. NaN pixels (masked or UNSEEN
    regions of the map) are ignored.

    The function counts the number of times each pixel contributes (its number
    of valid, finite values across patches) and divides the sum by that count,
    rather than by the total number of patches. Pixels with fewer than
    ``min_count`` contributions are set to NaN.

    Parameters
    ----------
    patches : sequence of numpy.ndarray
        Patches of identical shape, e.g. the output of :func:`extract_patches`.
    min_count : int or float, optional
        Minimum number of contributions for a pixel to be kept. A float in
        (0, 1) is a fraction of the number of patches. Default 1 (keep every
        pixel with at least one contribution).
    return_counts : bool, optional
        If True, also return the per-pixel count map. Default False.

    Returns
    -------
    stacked : numpy.ndarray
        The mean 2D patch. Display with ``plt.imshow(stacked)``.
    counts : numpy.ndarray of int
        Number of patches contributing to each pixel. Only returned if
        ``return_counts`` is True.
    """
    stack = np.array(patches, dtype=float)
    if isinstance(min_count, float) and 0.0 < min_count < 1.0:
        min_count = int(np.ceil(min_count * len(stack)))
    min_count = max(int(min_count), 1)

    valid = np.isfinite(stack)
    counts = valid.sum(axis=0)
    total = np.where(valid, stack, 0.0).sum(axis=0)

    mean_stacked = np.full(counts.shape, np.nan)
    np.divide(total, counts, out=mean_stacked, where=counts >= min_count)

    if return_counts:
        return mean_stacked, counts
    return mean_stacked


# ---------------------------------------------------------------------------
# Profile characterisation
# ---------------------------------------------------------------------------
def radial_profile(stacked, reso_arcmin=3.0, n_bins=None):
    """Azimuthally average a stacked patch into a 1D radial profile.

    Collapses the 2D stacked image to value-versus-radius by averaging in
    concentric annuli about the centre. This 1D profile is the characterisation
    of the mean peak: a central maximum, and (for CMB temperature) a faint
    acoustic ring further out.

    Parameters
    ----------
    stacked : numpy.ndarray
        Square 2D stacked patch from :func:`stack_patches`.
    reso_arcmin : float, optional
        Pixel size in arcminutes, so the returned radius is in physical angular
        units. Default 3.0.
    n_bins : int, optional
        Number of radial bins. Default is half the patch side length.

    Returns
    -------
    radius_arcmin : numpy.ndarray
        Bin-centre radii in arcminutes.
    profile : numpy.ndarray
        Mean value in each annulus.
    """
    ny, nx = stacked.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0

    y, x = np.indices(stacked.shape)
    r_pix = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)

    if n_bins is None:
        n_bins = nx // 2

    r_max = nx // 2
    bin_edges = np.linspace(0, r_max, n_bins + 1)
    bin_centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    profile = np.full(n_bins, np.nan)
    flat_r = r_pix.ravel()
    flat_v = stacked.ravel()
    for i in range(n_bins):
        in_bin = (flat_r >= bin_edges[i]) & (flat_r < bin_edges[i + 1])
        if np.any(in_bin):
            profile[i] = np.nanmean(flat_v[in_bin])

    radius_arcmin = bin_centres * reso_arcmin
    return radius_arcmin, profile