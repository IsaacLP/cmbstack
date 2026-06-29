[![Made at Code/Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)](https://semaphorep.github.io/codeastro/)
[![PyPI version](https://img.shields.io/pypi/v/cmbstack.svg)](https://pypi.org/project/cmbstack/)
[![Python version](https://img.shields.io/pypi/pyversions/cmbstack.svg)](https://pypi.org/project/cmbstack/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Development Status](https://img.shields.io/badge/Status-In_Development-brightgreen.svg)

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