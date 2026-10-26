import numpy as np
import healpy as hp

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

    lon_expected, lat_expected = hp.pix2ang(nside, spike_pix, lonlat=True)
    lon_found, lat_found = positions[0]

    assert np.isclose(lon_found, lon_expected)
    assert np.isclose(lat_found, lat_expected)


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


def test_peaks_near_mask_dropped():
    '''Test if find_peaks discards peaks touching the mask, or within mask_distance_deg of it'''
    nside = 64
    m = np.zeros(hp.nside2npix(nside))

    edge_pix = hp.ang2pix(nside, 0.0, 0.0, lonlat=True)
    far_pix = hp.ang2pix(nside, 180.0, 0.0, lonlat=True)
    mid_pix = hp.ang2pix(nside, 0.0, 10.0, lonlat=True)
    m[[edge_pix, far_pix, mid_pix]] = 10.0

    # Mask a neighbour of edge_pix, plus a pixel 5 deg away from mid_pix.
    m[hp.get_all_neighbours(nside, edge_pix)[0]] = hp.UNSEEN
    m[hp.ang2pix(nside, 0.0, 15.0, lonlat=True)] = hp.UNSEEN

    positions, n_rejected = find_peaks(m, nside=nside, threshold=5.0, return_n_rejected=True)
    found = hp.ang2pix(nside, positions[:, 0], positions[:, 1], lonlat=True)
    assert set(found) == {far_pix, mid_pix}
    assert n_rejected == 1

    positions, n_rejected = find_peaks(m, nside=nside, threshold=5.0, mask_distance_deg=6.0,
                                       return_n_rejected=True)
    found = hp.ang2pix(nside, positions[:, 0], positions[:, 1], lonlat=True)
    assert set(found) == {far_pix}
    assert n_rejected == 2


if __name__ == '__main__':

    test_high_threshold()

    test_high_n_peaks()

    test_single_peak()