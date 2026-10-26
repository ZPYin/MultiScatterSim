# MUSS–MSCART Simulation Comparison

This directory compares MUSS simulation results with reference data generated using the third-party MSCART model for ground-based and spaceborne lidar cases. The comparison covers parallel and perpendicular polarization signals, linear depolarization ratios for individual range bins, and cumulative linear depolarization ratios.

MSCART is a third-party model. Its source code and executables are not distributed with this repository. The MSCART results provided here serve as reference data for evaluating MUSS under the corresponding case configurations.

This directory contains saved comparison data, a plotting script, dependency requirements, and data provenance records.

The two cases are organized as follows:

| Directory | Comparison data |
| --- | --- |
| [ground](ground/) | MUSS results and MSCART reference data for ground-based lidar under low and high extinction conditions, with full fields of view of 0.5 and 2.0 mrad. |
| [satellite](satellite/) | MUSS results and MSCART reference data for the Hu2001 spaceborne lidar case, with full fields of view of 0.001, 0.13, and 13 mrad. |

See the [ground-based comparison](ground/comparison.png), the [spaceborne comparison](satellite/comparison.png), and `summary.csv` in each case directory. PDF versions are available for printing. Each `comparison.csv` contains the mean signals in both polarization channels, the depolarization ratio for each range bin, and the cumulative depolarization ratio used in the figures.

## Signal Definitions and Comparison Settings

The linear polarization channels are defined as `I∥=(I+Q)/2` and `I⊥=(I−Q)/2`, with the linear depolarization ratio `δ=I⊥/I∥`. Each channel is averaged across independent repeats before accumulation and ratio calculation. Shaded bands show the range across independent repeats.

- **Ground-based case:** 355 nm wavelength and 5 m range bins. Summary statistics use raw signals over cloud penetration depths of 20–120 m. MSCART outputs per km² are divided by 10⁶ to convert them to outputs per m². The figures show cumulative results up to 150 m.
- **Spaceborne case:** 532 nm wavelength and 3 m range bins. Summary statistics use range-corrected signals over cloud penetration depths of 0–150 m and 0–300 m. The figures show cumulative results up to 300 m.
- **Photon counts:** Each MUSS repeat uses 2 million photons. Each MSCART repeat uses 2 million photons for the ground-based case and 10 million photons for the spaceborne case. Contributions from scattering orders 4 and higher are obtained by subtracting the first three orders from the total signal.

## Plot Comparisons from Saved Data

Both the MUSS simulation results and the MSCART reference data are included. The plotting script reads these saved data directly to generate comparison figures.

In PyCharm, open `compare.py` in this directory, right-click, and select **Run**. Alternatively, execute `python compare.py` from this directory. By default, the script reads `model_results.csv` and the MSCART reference result files for both cases and generates the ground-based and spaceborne comparison figures. It prints the input CSV paths and output figure paths to the console. On Windows, the PNG images open automatically after plotting.

Use `--case ground` or `--case satellite` to plot a single case. Add `--no-show` to save the output files without opening the images.

Run the following commands from this directory:

```bash
python -m pip install -r requirements.txt
python compare.py --case ground
python compare.py --case satellite
```

Plotting requires only Python and the dependencies listed in `requirements.txt`. MATLAB and MSCART are not required.

## Reference Data Sources

The MSCART reference results are stored in `ground/low_result.nc` and `ground/high_result.nc` for the ground-based case, and in `satellite/result.nc` for the spaceborne case. These files contain results generated using MSCART for comparison with MUSS simulation results.

`source_info.json` records the reference data sources, SHA256 checksums, and data verification details.
