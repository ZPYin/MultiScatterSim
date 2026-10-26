# MSMC Optical Preprocessing

This directory contains the cloud optical preprocessing shared by the CPU and GPU semi-analytical polarized Monte Carlo models. Before photon transport, it converts cloud microphysical parameters and a Mie scattering database into a read-only `optics` structure that both backends can use directly.

## Module Scope

```text
Mie scattering database
+ lidarCfg.lambda
+ cloudCfg.clReff, N_d, kappa
          │
          ▼
prepareMSMCOptics
          │
          ├──► MSMCSim       (CPU photon transport)
          └──► MSMCSimGPU    (GPU photon transport)
```

This module does not perform photon free propagation, multiple scattering loops, Stokes reference-plane rotations, return-signal accumulation by range bin, or microphysical parameter retrieval.

The repository includes a scattering database for the minimal examples. The default file is `data/water_cloud_scattering_phase_function_example.mat`. Other compatible high-resolution databases can be specified explicitly using `scatMatrixFile`.

## Files

### `prepareMSMCOptics.m`

Calls `cloudScatteringMatrix` to calculate bulk optical properties for the specified wavelength and droplet size distribution, producing:

- Layer-specific extinction coefficients `muExt`;
- Layer-specific scattering coefficients `muSca`;
- Single-scattering albedo `albedo`;
- Mueller matrix elements `S11`, `S12`, `S33`, and `S34`;
- Scattering-angle grid `theta`;
- `cdfS11` and `cdfS12` for scattering-angle sampling.

### `prepareMSMCOpticsCached.m`

Provides memory and disk caching for batch forward simulations and repeated runs. For a fixed wavelength, `Reff/kappa` profiles, and normalized `N_d(z)` profile shape, changing the absolute `N_d` scale changes only `muExt` and `muSca` linearly. The Mueller matrices and normalized scattering-angle distributions can be reused.

Typical uses include:

- Large `Reff-N_d` parameter sweeps;
- Repeated Monte Carlo simulations with different random seeds;
- Batch calculations at multiple cloud base heights using the same in-cloud microphysical profiles.

## Usage

Use the project setup script to load all shared dependencies:

```matlab
setupMSToolbox
```

Alternatively, recursively add the entire `lib` directory:

```matlab
addpath(genpath(fullfile(projectRoot, 'lib')));
```

Standard calls:

```matlab
optics = prepareMSMCOptics(lidarCfg, cloudCfg, ...
    'scatMatrixFile', scatteringDatabase);

[SVcpu, diagCPU] = MSMCSim(lidarCfg, cloudCfg, ...
    'preparedOptics', optics);
[SVgpu, diagGPU] = MSMCSimGPU(lidarCfg, cloudCfg, ...
    'preparedOptics', optics);
```

Cached call:

```matlab
[optics, cacheInfo] = prepareMSMCOpticsCached( ...
    lidarCfg, cloudCfg, ...
    'scatMatrixFile', scatteringDatabase);
```

The CPU and GPU backends must share the same optical definitions. Calculations of extinction coefficients, Mueller matrices, and scattering-angle CDFs must not be maintained separately in the two backends.
