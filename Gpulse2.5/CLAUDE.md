# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Gpulse is a C++ library for calculating tightly focused laser pulse fields using Stratton-Chu diffraction integrals. It models radially polarized light focused by a parabolic mirror with ultra-short pulse temporal characteristics. The library uses memory-resident pre-computation for fast field queries at arbitrary space-time coordinates.

## Build and Run Commands

### Building
```bash
mkdir build && cd build
cmake ..                    # Configure with CMake (requires CMake 3.10+)
make                        # Build library and tools
```

The build creates:
- `libgpulse_lib.so`: Shared library
- `test_gpulse`: Basic field query test
- `export_field`: HDF5 field snapshot exporter (requires HDF5)
- `export_temporal`: Temporal intensity profiler (requires HDF5)
- `plot_field.py`, `plot_temporal.py`: Python visualization scripts

### Running Tests
```bash
cd build
./test_gpulse              # Run basic field calculation test
```

### Running Export Tools
```bash
cd build
./export_field [t_seconds]     # Export 2D field snapshot at time t (default t=0)
./export_temporal               # Export temporal intensity profile (-50fs to +50fs)
```

### Dependencies
- **Required**: OpenMP (for parallel processing)
- **Optional**: HDF5 C++ (for export tools only)
- **C++ Standard**: C++17
- **Compiler flags**: `-O3 -march=native -ffast-math`

## Code Architecture

### Core Components

#### 1. Field Solver (`src/Gpulse.cpp`)
The solver implements the Stratton-Chu vector diffraction integral for electromagnetic fields:

- **Pre-computation phase**: `run_solver()` calculates field values across the focal plane at multiple frequencies
- **Memory caching**: Results stored in `.gpulse_cache.bin` to avoid recomputation (validates parameters match)
- **Optimization**: Pre-calculates source-dependent terms (N×B, N·E) once per frequency before parallel pixel loops

#### 2. SpectralField Class (internal)
Stores pre-calculated field data in a 3D grid structure:
- **Dimensions**: `[pixel_x, pixel_y, frequency]`
- **Data**: Complex values for 6 field components (Ex, Ey, Ez, Bx, By, Bz)
- **Interpolation**: Bilinear interpolation for continuous (x,y) queries

#### 3. Public API (`include/Gpulse.h`)
- `Gpulse_Init()`: Initialize solver (lazy initialization, called automatically on first query)
- `Gpulse_Free()`: Release pre-computed memory
- `Gpulse(t, x, y, z, ...)`: Query real-valued fields
- `GpulseComplex(t, x, y, z, ...)`: Query complex-valued fields (for projection calculations)

#### 4. 2.5D Approximation
The z-dependence is separated using propagation phase `exp(i(kz - ωt))`. The transverse field distribution is assumed constant within the Rayleigh range. This is accurate for the focal depth region.

### Physical Parameters (`include/constants.h`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| LAMBDA0 | 1.8 μm | Center wavelength |
| W0 | 8.9 mm | Gaussian beam waist |
| P_PEAK | 0.3 TW | Peak laser power |
| WIDTH | 700 nm | Spectral bandwidth (FWHM) |
| ORDER | 7 | Super-Gaussian order |
| F0 | 6.35 mm | Parabolic mirror focal length |
| D | 25.4 mm | Mirror diameter |
| HOLE | 5 mm | Central obscuration diameter |
| RANGE | 18 μm | Observation window half-width |
| N_FREQ | 101 | Frequency sampling points |
| local_ycells/zcells | 100 | Grid resolution |

### Polarization Modes

The `RADIAL` constant (set to 1) enables Mosaic phase plate approximation for radial polarization:
- **Linear mode** (`RADIAL=0`): Standard linear polarization
- **Radial mode** (`RADIAL=1`): 4-quadrant phase plate approximation (π/4 rotation per quadrant)

### Numerical Integration

- **Source grid**: Simpson's rule on parabolic mirror surface
- **Spectral integration**: Simpson's rule over frequency domain
- **Resolution**: `NI=8` → (2^8+1) radial points, (2^9+1) angular points

## Key Implementation Details

### Cache System
The solver checks for `.gpulse_cache.bin` before computing. The cache header contains all physical parameters—if any parameter changes, the cache is invalidated and recomputed.

### Thread Safety
- `Gpulse_Init()` uses `std::call_once` for safe lazy initialization
- Main field computation loop uses OpenMP parallel for with dynamic scheduling

### Boundary Handling
- Queries outside `[±RANGE, ±RANGE]` return zero fields
- Bilinear interpolation handles off-grid points by clamping to valid range

### Field Reconstruction
Time-domain fields are reconstructed via inverse Fourier transform:
```
E(t) = Σ_n E_n(ω_n) × exp(-i×ω_n×t)
```
The phase factor includes both temporal oscillation and z-propagation:
```
phase = ω_n×(t - z/c)
```

## Output Data

### field_snapshot.h5 (export_field)
Datasets: `Ex_r, Ex_i, Ey_r, Ey_i, Ez_r, Ez_i, Bx_r, Bx_i, By_r, By_i, Bz_r, Bz_i`
Attributes: `range` (spatial extent), `t` (time), `radial` (polarization mode)

### temporal_intensity.h5 (export_temporal)
Datasets: `time` (seconds), `intensity` (W/m²)
Attribute: `radial` (polarization mode)

## Important Notes

1. **Coordinate system**: z=0 is the focal plane. Positive z is propagation direction.
2. **Time reference**: t=0 corresponds to pulse peak arrival at focal plane.
3. **Units**: All values in SI units (meters, seconds, Tesla, V/m).
4. **Memory**: Pre-computed field data occupies significant memory (~hundreds of MB depending on resolution).
5. **Performance**: First initialization takes seconds; subsequent queries are fast (microseconds).
