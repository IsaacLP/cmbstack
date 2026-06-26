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

# Function 1
def load_cl(path):
    """Load a power-spectrum file and return C_ell for the requested spectrum.

    The expected file columns are ell, Dl_TT, Dl_TE, Dl_EE, Dl_BB, Dl_dd, with
    D_ell in uK^2. The chosen column is converted from D_ell to C_ell via
    :func:`dl_to_cl` before being returned.

    Parameters
    ----------
    path : str
        Path to the whitespace-delimited spectrum file.
    column : str, optional
        Which spectrum to return: one of "TT", "TE", "EE", "BB", "dd".
        Default "TT".

    Returns
    -------
    cl : numpy.ndarray
        C_ell array indexed from ell=0, suitable for passing to
        :func:`simulate_map`.
    """
    ell, dl = np.loadtxt(path,usecols=(0,1),unpack=True)
    ell = ell.astype(int)
    lmax = ell.max()  
    cl = dl_to_cl(ell,dl,lmax)

    return cl

# Function 2
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


def normalize_map(m, remove_monopole=True):
    """Subtract the monopole and divide by the standard deviation.

    After this, peak thresholds can be expressed in units of sigma, which is the
    natural convention for peak statistics.

    Parameters
    ----------
    m : array_like
        Input HEALPix map. May contain UNSEEN/NaN pixels, which are ignored in
        the mean and standard deviation.
    remove_monopole : bool, optional
        If True, subtract the mean (monopole) before scaling. Default True.

    Returns
    -------
    m_norm : numpy.ndarray
        The normalized map, with (if monopole removed) mean ~0 and std ~1.
    """
    m_norm = hp.remove_monopole(m) / m.std()
    return m_norm


def load_map(path, field=0):
    """Wraps ``hp.read_map``"""
    return hp.read_map(path, field=field)
 
 
def save_map(path, sky_map, overwrite=True):
    """Wraps ``hp.write_map``"""
    hp.write_map(path, sky_map, overwrite=overwrite)