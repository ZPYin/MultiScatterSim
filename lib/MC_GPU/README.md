# GPU Monte Carlo Backend

This directory contains the GPU implementation of the semi-analytical polarized Monte Carlo model. The CPU and GPU backends share cloud optical properties, Mueller matrices, and physical definitions; both implement the same physical model.

## Public Entry Points

- CPU: `lib/MC/MSMCSim.m`
- GPU: `lib/MC_GPU/MSMCSimGPU.m`
- The GPU entry point forwards its arguments to `MSMCSimGPUBackend.m` in this directory.
- If a GPU is unavailable or a configuration is unsupported, `fallbackToCPU` controls whether execution falls back to `MSMCSim`.

The project setup script, `setupMSToolbox.m`, recursively adds `lib` to the MATLAB search path.

If you do not use the setup script, recursively add the entire `lib` directory. The GPU backend requires the shared optical preprocessing in `lib/MC_Optics` and the CPU fallback functions in `lib/MC`:

```matlab
addpath(genpath(fullfile(projectRoot, 'lib')));
```

Here, `projectRoot` is the path to the repository root.

## Corresponding CPU and GPU Physical Processes

| CPU implementation | GPU implementation | Physical process |
| --- | --- | --- |
| `freeRunning.m` | `freeRunningGPU.m` | Random free propagation in a layered medium |
| `rotSVPhi.m` | `rotSVPhiGPU.m` | Rotation of the Stokes reference plane |
| `scaRotAng.m` | `scaRotAngGPU.m` | Scattering angle and the two reference-plane rotation angles |
| `updateDir.m` | `updateDirGPU.m` | Propagation direction update after scattering |
| `sampleScatteringAngle.m` | `sampleScatteringAngleGPU.m` | Sampling from the scattering-angle distribution |

GPU functions process batches of photon states stored as `gpuArray` objects. Their CPU counterparts process individual photons, so their calling interfaces and array dimensions require adaptation for batched GPU execution.

## Internal Structure

```text
lib/
├── MC/                          CPU photon transport
├── MC_Optics/                   Shared CPU/GPU optical preprocessing
└── MC_GPU/
    ├── MSMCSimGPU.m             Public GPU entry point
    ├── MSMCSimGPUBackend.m      GPU scheduling and result aggregation
    ├── simulateLayeredBatchGPU.m
    ├── freeRunningGPU.m
    ├── rotSVPhiGPU.m
    ├── scaRotAngGPU.m
    ├── updateDirGPU.m
    ├── sampleScatteringAngleGPU.m
    └── README.md
```

- `MSMCSimGPUBackend.m` handles argument parsing, device selection, optical data transfer to the GPU, batch scheduling, result aggregation, diagnostic output, and CPU fallback.
- `simulateLayeredBatchGPU.m` executes the complete photon scattering loop, including free propagation, scattering sampling, Mueller matrix operations, receiver checks, and signal accumulation.

## Photon State Layouts

`simulateLayeredBatchGPU.m` supports two photon state layouts, selected with the `stateLayout` parameter:

- **`legacy` (default):** Keeps fixed-length arrays of current photon states for the entire batch and uses an active-photon mask to select photons for further propagation. Fixed array positions make it easier to compare individual photons with the original implementation and debug state updates.
- **`packed`:** Removes the states of photons that have left the cloud, been absorbed, or otherwise terminated, retaining only active photon states at each step. This layout aims to reduce GPU indexing writebacks and CPU–GPU synchronization overhead. Its performance depends on the GPU, batch size, and diagnostic settings.

Both layouts use the same physical model and can be used for simulations. Diagnostic collection is controlled separately by `diagnosticsLevel`. Neither layout automatically records complete photon trajectory histories or exports per-photon state logs.

The GPU backend requires MATLAB Parallel Computing Toolbox and a GPU supported by MATLAB.
