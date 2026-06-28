[![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)](https://semaphorep.github.io/codeastro/)

# cmbstack

Stack patches of the Cosmic Microwave Background (CMB) temperature sky. Go [here](https://cmbstack-docs.readthedocs.io/en/latest/index.html) to see the full documentation and a worked example.

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
from cmbstack.main import StackingPipeline

pipeline = StackingPipeline.from_cl("path/to/spectrum.cl", nside=1024, seed=42)
StackedResult = pipeline.run()
```

**From an existing FITS map:**
```python
pipeline = StackingPipeline.from_fits("path/to/map.fits", field=0)
StackedResult = pipeline.run()
```

**From a map array already in memory:**
```python
pipeline = StackingPipeline.from_map(sky_map)
StackedResult = pipeline.run()
```