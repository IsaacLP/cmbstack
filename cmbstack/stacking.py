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
def find_peaks(sky_map, nside, threshold=None, n_peaks=None):
    """Find local maxima of a HEALPix map and return their sky positions.

    A pixel is a local maximum if its value is strictly greater than all of its
    immediate HEALPix neighbours. Peaks can be filtered by a significance
    threshold and/or capped at the ``n_peaks`` highest.

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
        threshold). Default None (keep all).

    Returns
    -------
    positions : numpy.ndarray, shape (N, 2)
        Sky positions of the selected peaks as ``(theta, phi)`` in radians, the
        same format accepted by :func:`extract_patches`.
    """
    # hp.hotspots returns (max_map, minima_pix, maxima_pix); we want the maxima.
    _, _, maxima_pix = hp.hotspots(sky_map)
    maxima_pix = np.asarray(maxima_pix)

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

    theta, phi = hp.pix2ang(nside, maxima_pix)
    return np.column_stack([theta, phi])


# ---------------------------------------------------------------------------
# Patch extraction
# ---------------------------------------------------------------------------
def extract_patches(sky_map, positions, size_deg=10.0, reso_arcmin=3.0):
    """Extract fixed-grid gnomonic patches centred on each position.

    Each patch is a square 2D array produced by a gnomonic (tangent-plane)
    projection centred on the position, so every patch shares the same grid and
    the centre pixel always corresponds to the position itself. 

    Parameters
    ----------
    sky_map : numpy.ndarray
        Input HEALPix map (RING ordering). Any scalar field.
    positions : array_like, shape (N, 2)
        Sky positions as ``(theta, phi)`` in radians (e.g. the output of
        :func:`find_peaks`, or an external catalogue converted to this format).
    size_deg : float, optional
        Full side length of the square patch in degrees. Default 10.0.
    reso_arcmin : float, optional
        Pixel size of the projected patch in arcminutes. Default 3.0.

    Returns
    -------
    patches : list of numpy.ndarray
        One square 2D array per position, all of identical shape
        ``(xsize, xsize)`` with ``xsize = size_deg * 60 / reso_arcmin``.
    """
    positions = np.atleast_2d(positions)
    xsize = int(round(size_deg * 60.0 / reso_arcmin))

    patches = []
    for theta, phi in positions:
        # gnomview's rot expects (lon, lat) in degrees.
        lon = np.degrees(phi)
        lat = 90.0 - np.degrees(theta)
        patch = hp.gnomview(
            sky_map,
            rot=(lon, lat),
            xsize=xsize,
            reso=reso_arcmin,
            return_projected_map=True,
            no_plot=True,
        )
        patches.append(np.asarray(patch))

    return patches


# ---------------------------------------------------------------------------
# Stacking
# ---------------------------------------------------------------------------
def stack_patches(patches):
    """Average patches pixel-by-pixel into a single stacked image.

    Because every patch shares the same fixed grid (see :func:`extract_patches`),
    the mean is a genuine stacked image: incoherent noise averages towards zero
    while the coherent central profile survives. NaN/UNSEEN pixels (patch edges
    that fall off the map near the poles) are ignored.

    Parameters
    ----------
    patches : sequence of numpy.ndarray
        Patches of identical shape, e.g. the output of :func:`extract_patches`.

    Returns
    -------
    stacked : numpy.ndarray
        The mean 2D patch. Display with ``plt.imshow(stacked)``.
    """
    stack = np.array(patches, dtype=float)
    mean_stacked = np.nanmean(stack, axis=0)
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