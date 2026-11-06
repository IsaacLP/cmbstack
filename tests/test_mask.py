import warnings

import numpy as np
import healpy as hp
import pytest

from cmbstack import StackingPipeline
from cmbstack.maps import normalize_map
from cmbstack.stacking import stack_patches, extract_patches, iter_patches


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


def test_run_rectangular_patches():
    '''Test if run() passes size_deg_y through to the patches and the stack'''
    _, m_masked = _masked_map()
    pipeline = StackingPipeline.from_map(m_masked).run(size_deg=10.0, size_deg_y=6.0,
                                                       reso_arcmin=30.0, threshold=2.0,
                                                       profile=False)

    # extract_patches returns (ysize, xsize) arrays.
    assert pipeline.patches[0].shape == (12, 20)
    assert pipeline.stacked.shape == (12, 20)


def test_stack_patches_counts():
    '''Test if stack_patches divides by the per-pixel count and applies min_count'''
    nan = np.nan
    patches = [
        np.array([[1.0, 2.0], [nan, nan]]),
        np.array([[3.0, nan], [5.0, nan]]),
        np.array([[5.0, 4.0], [nan, nan]]),
    ]

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no "Mean of empty slice" for the all-NaN pixel
        stacked, counts = stack_patches(patches, return_counts=True)

    assert np.array_equal(counts, [[3, 2], [1, 0]])
    assert np.allclose(stacked, [[3.0, 3.0], [5.0, nan]], equal_nan=True)

    # Pixels seen by fewer than 2 patches are dropped, also as a fraction.
    expected = [[3.0, 3.0], [nan, nan]]
    assert np.allclose(stack_patches(patches, min_count=2), expected, equal_nan=True)
    assert np.allclose(stack_patches(patches, min_count=0.5), expected, equal_nan=True)


def _masked_patch_args():
    '''Masked map, random positions and angles for the streaming tests'''
    _, m_masked = _masked_map()
    rng = np.random.default_rng(1)
    positions = np.column_stack([rng.uniform(0, 360, 50), rng.uniform(-60, 60, 50)])
    psi = rng.uniform(0, 180, 50)
    return m_masked, positions, dict(psi=psi, size_deg=10.0, size_deg_y=6.0, reso_arcmin=30.0)


def test_iter_patches_matches_extract_patches():
    '''Test if iter_patches yields the same patches as extract_patches, in the same order'''
    m_masked, positions, kwargs = _masked_patch_args()

    listed = extract_patches(m_masked, positions, **kwargs)
    streamed = list(iter_patches(m_masked, positions, **kwargs))

    assert len(streamed) == len(listed) == len(positions)
    for a, b in zip(listed, streamed):
        assert np.array_equal(a, b, equal_nan=True)


def test_stack_patches_generator_matches_list():
    '''Test if stacking a generator gives the same stack and counts as stacking a list'''
    m_masked, positions, kwargs = _masked_patch_args()
    patches = extract_patches(m_masked, positions, **kwargs)
    assert any(np.isnan(p).any() for p in patches)  # the mask reaches some patches

    for min_count in [1, 0.9]:
        s_list, c_list = stack_patches(patches, min_count=min_count, return_counts=True)
        s_gen, c_gen = stack_patches(iter_patches(m_masked, positions, **kwargs),
                                     min_count=min_count, return_counts=True)
        assert np.array_equal(c_list, c_gen)
        assert np.allclose(s_list, s_gen, equal_nan=True)


def test_stack_patches_errors():
    '''Test if stack_patches rejects an empty input and patches of different shapes'''
    with pytest.raises(ValueError):
        stack_patches(iter([]))
    with pytest.raises(ValueError):
        stack_patches([np.zeros((2, 2)), np.zeros((2, 3))])


def test_run_keep_patches_false():
    '''Test if run(keep_patches=False) gives the same stack and counts without storing the patches'''
    _, m_masked = _masked_map()
    kwargs = dict(size_deg=10.0, reso_arcmin=30.0, threshold=2.0)

    kept = StackingPipeline.from_map(m_masked).run(**kwargs)
    streamed = StackingPipeline.from_map(m_masked).run(keep_patches=False, **kwargs)

    assert streamed.patches is None
    assert len(kept.patches) == len(streamed.positions)
    assert np.array_equal(kept.counts, streamed.counts)
    assert np.allclose(kept.stacked, streamed.stacked, equal_nan=True)
    assert np.allclose(kept.profile, streamed.profile, equal_nan=True)


if __name__ == "__main__":
    pytest.main([__file__])
