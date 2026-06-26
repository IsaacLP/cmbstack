import sys
import numpy as np
import healpy as hp

sys.path.append("..")

from cmbstack.stacking import find_peaks
from cmbstack.maps import normalize_map


def test_single_peak():
    '''Test if find_peaks returns the corrected coordinates for a known peak'''
    nside = 64
    n_pix = hp.nside2npix(nside)
    m = np.zeros(n_pix)

    spike_pix = n_pix // 2
    m[spike_pix] = 10.0

    positions = find_peaks(m, nside=nside, threshold=5.0)

    assert len(positions) == 1

    theta_expected, phi_expected = hp.pix2ang(nside, spike_pix)
    theta_found, phi_found = positions[0]

    assert np.isclose(theta_found, theta_expected)
    assert np.isclose(phi_found, phi_expected)


def test_high_threshold():
    '''Test if find_peaks returns empty array for very high threshold'''
    cl = np.zeros(384)
    cl[2:] = 1.0/np.arange(2,384)**2
    np.random.seed(0)
    m = hp.synfast(cl, nside=128)
    m_norm = normalize_map(m)
    threshold = 1e3

    positions = find_peaks(m_norm,nside=128,threshold=threshold)

    assert positions.size == 0


def test_high_n_peaks():
    '''Test if find_peaks returns all the peaks when n_peaks is greater than the total number of found peaks'''
    cl = np.zeros(384)
    cl[2:] = 1.0/np.arange(2,384)**2
    np.random.seed(0)
    m = hp.synfast(cl, nside=128)

    positions = find_peaks(m,nside=128)

    n_peaks = len(positions)

    test_n_peaks = n_peaks + 1

    new_positions = find_peaks(m,nside=128,n_peaks=test_n_peaks)

    assert len(new_positions) == len(positions)


if __name__ == '__main__':

    test_high_threshold()

    test_high_n_peaks()

    test_single_peak()