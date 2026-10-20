"""Map simulation and preprocessing utilities for HEALPix CMB maps.

Provides the building blocks for turning a raw power spectrum into a
normalised HEALPix map ready for stacking:

  1. :func:`load_cl`      — read a D_ell spectrum file and convert to C_ell
  2. :func:`simulate_map` — draw a Gaussian random realisation with ``healpy.synfast``
  3. :func:`normalize_map`— subtract the monopole and divide by the std

These functions are intentionally field-agnostic: they work on any scalar
power spectrum (temperature TT, lensing convergence, y-map, ...).
"""

import numpy as np
import healpy as hp
import pandas as pd


def dl_to_cl(ell, dl,lmax):
    """Convert D_ell = ell(ell+1) C_ell / (2 pi) to C_ell.

    Parameters
    ----------
    ell : array_like
        Multipole values. May start at 0; the ell=0 and ell=1 entries are set
        to zero in the output to avoid division by zero (they carry no usable
        power for this purpose).
    dl : array_like
        D_ell values in the same units you want C_ell returned in (e.g. uK^2).

    Returns
    -------
    cl : numpy.ndarray
        The angular power spectrum C_ell, same shape as ``dl``.

    Notes
    -----
    The inverse normalization is ``C_ell = D_ell * 2*pi / (ell*(ell+1))``.
    """
    # Convert D_ell -> C_ell
    cl_vals = 2.0 * np.pi * dl / (ell * (ell + 1.0))

    # Build a full array indexed from ell=0, with 0,1 set to zero
    lmax = np.max(ell)
    cl = np.zeros(lmax + 1)
    cl[ell] = cl_vals

    return cl


def load_cl(path, column="TT"):
    """Load a power-spectrum file and return C_ell for the requested spectrum.

    The expected file columns are ell, Dl_TT, Dl_TE, Dl_EE, Dl_BB, Dl_pp, with
    D_ell in uK^2 (with names L, TT, TE, EE, BB, PP).
    The chosen column is converted from D_ell to C_ell via
    :func:`dl_to_cl` before being returned.

    Parameters
    ----------
    path : str
        Path to the whitespace-delimited spectrum file.
    column : str, optional
        Which spectrum to return: one of "TT", "TE", "EE", "BB", "PP".
        Default "TT".

    Returns
    -------
    cl : numpy.ndarray
        C_ell array indexed from ell=0, suitable for passing to
        :func:`simulate_map`.
    """

    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["L", "TT", "TE", "EE", "BB", "PP"],
    )
    ell = df['L'].values.astype(int)
    dl = df[column].values
    lmax = ell.max()
    cl = dl_to_cl(ell, dl, lmax)

    return cl


def simulate_map(cl, nside=128, seed=None):
    """Simulate a Gaussian random HEALPix map from a power spectrum.

    Thin wrapper over ``healpy.synfast`` with an optional seed so results are
    reproducible in tests.

    Parameters
    ----------
    cl : array_like
        Angular power spectrum C_ell (not D_ell).
    nside : int, optional
        HEALPix resolution parameter. Default 128.
    seed : int or None, optional
        Seed for the random number generator. If None, the draw is random.

    Returns
    -------
    m : numpy.ndarray
        A HEALPix map (RING ordering) of length ``12 * nside**2``.
    """
    if seed:
        np.random.seed(seed)

    map = hp.synfast(cl, nside=nside) 

    return map


def normalize_map(m):
    """Subtract the monopole and divide by the standard deviation.

    After this, peak thresholds can be expressed in units of sigma, which is the
    natural convention for peak statistics.

    Parameters
    ----------
    m : array_like
        Input HEALPix map. May contain UNSEEN/NaN pixels or be a masked array
        (e.g. from ``hp.ma``); these pixels are ignored in the mean and
        standard deviation and set to UNSEEN in the output.

    Returns
    -------
    m_norm : numpy.ndarray
        The normalized map, with mean ~0 and std ~1.
    """
    m = np.ma.filled(np.ma.asarray(m, dtype=float), hp.UNSEEN)
    # mask_bad compares to UNSEEN with a tolerance, so UNSEEN stored as float32
    # (e.g. Planck maps from hp.read_map) is still recognised; it also flags NaN/inf.
    good = ~hp.mask_bad(m)
    m_norm = np.full_like(m, hp.UNSEEN)
    m_good = m[good] - m[good].mean()
    m_norm[good] = m_good / m_good.std()
    return m_norm


def load_map(path, field=0):
    """Wraps ``hp.read_map``"""
    return hp.read_map(path, field=field)
 
 
def save_map(path, sky_map, overwrite=True):
    """Wraps ``hp.write_map``"""
    hp.write_map(path, sky_map, overwrite=overwrite)