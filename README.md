# CMB Peak Stacking Pipeline

## Overview

`cmbstack` is a Python package for stacking patches of the Cosmic Microwave Background (CMB) temperature sky. It accepts input as a theoretical power spectrum, a [HEALPix](https://healpy.readthedocs.io/en/latest/index.html) FITS file, or a map array already in memory. From there it detects local maxima, extracts gnomonic (flat-sky) patches around each peak, and averages them. This stacking procedure enhances the coherent peak profile while suppressing uncorrelated noise.

## Installation

```bash
pip install -e .
```

## Quick Start

**From a power spectrum file:**
```python
from cmbstack.main import StackingPipeline

pipeline = StackingPipeline.from_cl("path/to/spectrum.dat", nside=1024, seed=42)
pipeline.run()
```

**From an existing FITS map:**
```python
pipeline = StackingPipeline.from_fits("path/to/map.fits", field=0)
pipeline.run()
```

**From a map array already in memory:**
```python
pipeline = StackingPipeline.from_map(sky_map)
pipeline.run()
```

See `examples/from_theoretical_cl.ipynb` for a full worked example.

## Package Structure

```
cmbstack/
├── maps.py      — power-spectrum loading, map simulation, normalisation
├── stacking.py  — peak finding, patch extraction, stacking, radial profile
└── main.py      — StackingPipeline high-level class
```

## Workflow

> **Note:** Steps 1–2 apply when starting from a power spectrum. Use `StackingPipeline.from_fits` or `from_map` to skip them when working with a real or pre-simulated map.

### 1. Load the Power Spectrum — `maps.load_cl`

The input file contains the power spectrum as $D_\ell^{TT}$:

$$D_\ell \equiv \frac{\ell(\ell+1)}{2\pi} C_\ell$$

`load_cl` reads columns $(\ell, D_\ell)$ and converts to $C_\ell$:

$$C_\ell = \frac{2\pi}{\ell(\ell+1)} D_\ell, \qquad C_0 = C_1 = 0$$

---

### 2. Simulate a Sky Map — `maps.simulate_map`

A Gaussian random realization is drawn by sampling spherical harmonic coefficients $a_{\ell m}$ with variance $C_\ell$:

$$T(\hat{n}) = \sum_{\ell,m} a_{\ell m} \, Y_{\ell m}(\hat{n}), \qquad \langle |a_{\ell m}|^2 \rangle = C_\ell$$

This calls `healpy.synfast` internally. An optional `seed` makes runs reproducible.

---

### 3. Normalise the Map — `maps.normalize_map`

Before peak detection, the map is standardised so that thresholds have a clear statistical meaning:

$$T_{\text{norm}}(\hat{n}) = \frac{T(\hat{n}) - \langle T \rangle}{\sigma}$$

After this step the map has mean $\approx 0$ and standard deviation $= 1$, so peaks are measured in units of $\sigma$.

---

### 4. Detect Peaks — `stacking.find_peaks`

Local maxima are identified with `healpy.hotspots`: a pixel is a maximum if its value exceeds every immediate HEALPix neighbour. Peaks are filtered by a significance threshold $\nu$ (default $\nu = 3\sigma$) and optionally capped at the $N$ highest:

$$\text{Peaks} = \{\hat{n}_p \in \text{Maxima} \mid T_{\text{norm}}(\hat{n}_p) > \nu\}$$

Returns sky positions as $(θ, φ)$ pairs in radians.

---

### 5. Extract Patches — `stacking.extract_patches`

For each peak a square patch is cut using a gnomonic (tangent-plane) projection centred on $\hat{n}_p$. Every patch shares the same fixed pixel grid (side length `size_deg`, pixel scale `reso_arcmin`), so the centre pixel always corresponds to the peak itself and patches can be co-added directly.

---

### 6. Stack — `stacking.stack_patches`

Patches are averaged pixel-by-pixel:

$$S = \frac{1}{N} \sum_{i=1}^{N} P_i$$

Incoherent noise averages towards zero; the coherent central profile survives.

---

### 7. Radial Profile — `stacking.radial_profile`

The 2D stacked image is collapsed to a 1D profile by azimuthal averaging in concentric annuli about the centre. Returns bin-centre radii in arcminutes and the mean temperature in each annulus.

---

## Pipeline Constructors

`StackingPipeline` provides three entry points depending on where your data comes from:

| Constructor | Input | Notes |
|---|---|---|
| `from_cl(path, nside, seed)` | Power-spectrum file | Simulates a Gaussian random map via `healpy.synfast` |
| `from_fits(path, field=0)` | HEALPix FITS file | Loads the map with `maps.load_map`; `nside` is inferred automatically |
| `from_map(sky_map)` | NumPy array | Accepts any in-memory HEALPix map; `nside` is inferred automatically |

All three store the map in `pipeline.map` and share the same `run()` interface.

---

## Map I/O Utilities

`maps.load_map` and `maps.save_map` wrap the healpy FITS readers for convenience:

```python
from cmbstack import maps

m = maps.load_map("map.fits", field=0)   # wraps hp.read_map
maps.save_map("out.fits", m)             # wraps hp.write_map (overwrite=True by default)
```

---

## Pipeline Object

`StackingPipeline` stores every intermediate product as an attribute:

| Attribute | Content |
|---|---|
| `pipeline.map` | Raw simulated map |
| `pipeline.normalized` | Normalised map (units of $\sigma$) |
| `pipeline.positions` | Peak positions $(θ, φ)$ in radians |
| `pipeline.patches` | List of 2D gnomonic patches |
| `pipeline.stacked` | Mean stacked 2D image |
| `pipeline.radius` | Radial bin centres (arcmin) |
| `pipeline.profile` | Mean temperature per radial bin |
