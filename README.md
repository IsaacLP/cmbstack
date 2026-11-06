[![Made at Code/Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet)](https://semaphorep.github.io/codeastro/)
[![PyPI version](https://img.shields.io/pypi/v/cmbstack)](https://pypi.org/project/cmbstack/)
[![Python version](https://img.shields.io/pypi/pyversions/cmbstack)](https://pypi.org/project/cmbstack/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](https://opensource.org/licenses/MIT)
![Development Status](https://img.shields.io/badge/Status-In_Development-brightgreen)

# cmbstack

Stack patches of HEALPix maps, such as the CMB temperature, CMB lensing or galaxy density, around their peaks or around the positions of a catalogue. Go [here](https://cmbstack-docs.readthedocs.io/en/latest/index.html) to see the full documentation, and [here](https://github.com/IsaacLP/cmbstack/tree/main/docs/tutorials) for worked examples.

## Installation

Install the latest release from PyPI:

```bash
pip install cmbstack
```

Or install the development version from GitHub:

```bash
git clone https://github.com/IsaacLP/cmbstack.git
cd cmbstack
pip install -e .
```

## Quick Start

**From a power spectrum file:**
```python
from cmbstack import StackingPipeline

pipeline = StackingPipeline.from_cl("path/to/spectrum.cl", nside=1024, seed=42)
stacked_result = pipeline.run()
```

**From an existing FITS map:**
```python
pipeline = StackingPipeline.from_fits("path/to/map.fits", field=0)
stacked_result = pipeline.run()
```

**From a map array already in memory:**
```python
pipeline = StackingPipeline.from_map(sky_map)
stacked_result = pipeline.run()
```

A masked map (`healpy.ma`) keeps its mask: masked pixels are ignored, and `pipeline.counts` gives the number of patches contributing to each pixel.

**On the positions of a catalogue, with a rotation for each patch:**
```python
# positions: (N, 2) array of (lon, lat) in degrees; angles: (N,) rotations in degrees
pipeline = StackingPipeline.from_map(sky_map)
stacked_result = pipeline.run(positions=positions, psi=angles, size_deg=2.4, size_deg_y=12.0,
                              profile=False, keep_patches=False)
```

`keep_patches=False` stacks the patches one at a time instead of keeping them all in memory, which is needed for large catalogues.