import numpy as np
import healpy as hp


def dl_to_cl(ell, dl):
    # Convert D_ell -> C_ell
    cl_vals = 2.0 * np.pi * dl / (ell * (ell + 1.0))

    # Build a full array indexed from ell=0, with 0,1 set to zero
    lmax = np.max(ell)
    cl = np.zeros(lmax + 1)
    cl[ell] = cl_vals

    return cl


def load_cl(path):
    ell, dl = np.loadtxt(path,usecols=(0,1))
    cl = dl_to_cl(ell,dl)

    return cl


def simulate_map(cl, nside=128, seed=None):
    if seed:
        np.random.seed(seed)

    map = hp.synfast(cl, nside=nside) 

    return map


def normalize(map):
    return hp.remove_monopole(map) / map.std()