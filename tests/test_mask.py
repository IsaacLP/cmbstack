import numpy as np
import healpy as hp
import pytest

from cmbstack import StackingPipeline
from cmbstack.maps import normalize_map


def _masked_map(nside=32, dtype=np.float64):
    np.random.seed(0)
    m = np.random.normal(5.0, 2.0, hp.nside2npix(nside)).astype(dtype)
    m_masked = hp.ma(m)
    m_masked.mask = np.zeros(m.size, dtype=bool)
    m_masked.mask[: m.size // 4] = True
    return m, m_masked


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_normalize_map_keeps_mask(dtype):
    '''Test if normalize_map sets masked pixels to UNSEEN and ignores them in the statistics'''
    _, m_masked = _masked_map(dtype=dtype)
    mask = m_masked.mask

    m_norm = normalize_map(m_masked)

    assert np.all(m_norm[mask] == hp.UNSEEN)
    assert np.isclose(m_norm[~mask].mean(), 0.0, atol=1e-6)
    assert np.isclose(m_norm[~mask].std(), 1.0)


@pytest.mark.parametrize("dtype", [np.float64, np.float32])
def test_from_map_keeps_mask(dtype):
    '''Test if StackingPipeline.from_map preserves the mask of a masked array through run()'''
    _, m_masked = _masked_map(dtype=dtype)
    mask = m_masked.mask

    pipeline = StackingPipeline.from_map(m_masked)
    assert np.all(hp.mask_bad(pipeline.map) == mask)

    pipeline.run(threshold=2.0, profile=False)
    assert np.all(pipeline.normalized[mask] == hp.UNSEEN)

    # Peaks are found, and none of them inside the masked region.
    assert len(pipeline.positions) > 0
    pix = hp.ang2pix(pipeline.nside, pipeline.positions[:, 0], pipeline.positions[:, 1], lonlat=True)
    assert not np.any(mask[pix])

if __name__ == "__main__":
    pytest.main([__file__])
