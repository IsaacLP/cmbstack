import sys
import healpy as hp
import numpy as np
from pathlib import Path
import tempfile

sys.path.append("..")

from cmbstack.main import StackingPipeline
from cmbstack import maps

def test_fits_roundtrip(tmp_path):
    '''Test that the Stacking Pipeline gives the same result from a map saved in memory and loading a .fits file'''
    cl = np.zeros(384)
    cl[2:] = 1.0/np.arange(2,384)**2
    np.random.seed(0)
    m = hp.synfast(cl, nside=128)
    path = tmp_path / "m.fits"
    maps.save_map(str(path), m)

    direct = StackingPipeline.from_map(m).run(reso_arcmin=10, threshold=1.0, n_peaks=200)
    loaded = StackingPipeline.from_fits(str(path)).run(reso_arcmin=10, threshold=1.0, n_peaks=200)
    assert np.allclose(direct.stacked, loaded.stacked, equal_nan=True)

if __name__ == '__main__':
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        test_fits_roundtrip(tmp_path=tmp_path)