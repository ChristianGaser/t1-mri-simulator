# mri_simulate
Simulates T1-weighted MR images with optional atrophy, cortical thickness, WMHs, RF B1 inhomogeneities, and noise (Gaussian or Rician at a target WM SNR). Writes JSON sidecars with simulation metadata.

## Overview
`mri_simulate` generates simulated T1-weighted (T1w) images from a high-quality input (e.g., 0.5 mm Colin27 or a custom T1w). It can:

- Add noise: Gaussian (% of WM mean) or Rician at target WM SNR
- Apply RF B1 inhomogeneity (predefined A/B/C or simulated fields)
- Apply contrast change (power-law Y^x)
- Simulate white matter hyperintensities (WMHs)
- Simulate regional atrophy (atlas-based)
- Enforce a constant cortical thickness and produce a PVE-like segmentation used for synthesis

Thickness/PVE pipeline:
- GM is grown outward from WM using an Euclidean distance map to reach a specified thickness (global or 3-region, based on the Hammer atlas).
- To emulate partial volume effects, the label boundary is shifted over 20 offsets in [-0.25, 0.25] intensity range of the PVE label. Each offset yields a hard label (CSF=1, GM=2, WM=3); averaging across offsets produces a PVE-like label map.
- The resulting tissue maps replace SPM’s GM/WM/CSF posteriors in synthesis, using SPM’s Gaussian mixture parameters.

The function runs SPM12 segmentation automatically if needed and relies on CAT12 utilities. Ensure SPM12 and CAT12 are on your MATLAB path.

## Requirements
- MATLAB with SPM12 and CAT12 toolboxes in the path
- A T1-weighted NIfTI image (default examples use `colin27_t1_tal_hires.nii`)

## Inputs
### simu: Simulation parameters (struct)

Parameter | Description (Default)
----------|------------------------
name | Input image(s). A single T1w filename. (Default: `''`)
pn | Gaussian noise level as percent of the WM peak. (Default: `3`)
snrWM | If `>0`, add Rician magnitude noise with target SNR for white matter. Uses WM mean to derive complex noise sigma; when set, `pn` is ignored. (Default: `0`)
rng | RNG seed for reproducible noise; set `[]` for MATLAB default behavior. (Default: `0`)
contrast | Power-law exponent for contrast change. Image is normalized to [0,1], transformed as Y.^contrast, then rescaled to original min/max. Meaningful values to simulate contrast are 0.5 (low contrast) and 2 (high contrast). (Default: `1`)
resolution | Output voxel size: scalar (applied to x,y,z) or `[x y z]`. `NaN` keeps the original resolution. (Default: `NaN`)
WMH | Strength of white matter hyperintensities. `0`=off; `1`=mild; `2`=medium; `3`=strong; values `>=1` allowed. Larger values broaden the WMH prior via exponent `1/(WMH-0.8)` and scale the label contribution by `~1/WMH^0.75`. Constrained to (eroded) WM and modulated by a random field. (Default: `0`)
atrophy | Atrophy specification: `{atlasName, roiIds[], factors[]}`; factors >1 increase CSF (reduce GM) within ROIs. Either thickness or atrophy can be simulated. (Default: `[]`)
thickness | Cortical thickness in mm. Scalar = global; 3-vector = `[occipital rest frontal]` using Hammer atlas masks. Subcortical/cerebellar regions are excluded from thickness simulation and the original thickness values are kept. Either thickness or atrophy can be simulated. (Default: `0`)

### rf: RF bias field parameters (struct)

Parameter | Description (Default)
----------|------------------------
percent | Amplitude in percent; negative values invert the field. (Default: `20`)
type | `'A'|'B'|'C'` (predefined MNI fields) or numeric `[strength rngSeed]` for a simulated field. Strength in `1..4` (3–4 ~ stronger 7T-like). (Default: `[2 0]`)
save | Save the simulated bias field only when `type` is numeric; ignored for `'A'|'B'|'C'`. (Default: `0`)

## Defaults
If `simu` and/or `rf` are omitted or partially specified, missing fields are filled with defaults. If `simu.name` is empty, a file selection dialog opens.

```matlab
simu = struct('name', '', 'pn', 3, 'snrWM', 0, 'contrast', 1, ...
              'resolution', NaN, 'WMH', 0, 'atrophy', [], 'thickness', 0, 'rng', 0);
rf   = struct('percent', 20, 'type', [2 0], 'save', 0);
```

## Outputs
The function saves:
- Simulated image: `pn{pn}_{meanRes}mm_{name}{opts}.nii`
- Simulated masked image: `pn{pn}_{meanRes}mm_m{name}{opts}.nii`
- Ground-truth PVE label: `label_pve_{meanRes}mm_{name}{opts}.nii`
- If requested, RF field (simulated only): `{opts}_{meanRes}mm_{name}.nii`
- JSON sidecars (main and masked images): `{simuFile}.json` including tool metadata and SimulationParameters (voxel size, pn or snrWM, RF settings, thickness tag)

Notes:
- When `snrWM>0`, filenames use `snr{snrWM}` instead of `pn{pn}`.
- `{opts}` aggregates options, e.g., `_rf20_A`, `_WMH2`, `_hammers_28_2`, `_thickness1.5mm-2.5mm`.
- When thickness is used, the label is PVE-like from the boundary jittering averaging.
- When WMH is used, a 4th label contribution is added (WMH).

## Usage
```matlab
mri_simulate(simu, rf);
```

### 1b) Rician noise at target WM SNR
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'snrWM', 30, ...
              'pn', 0, 'resolution', NaN, 'rng', 0);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

## Examples

### 1) Basic simulation with specific noise and 0.5 mm voxels
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3, ...
              'resolution', 0.5, 'atrophy', [], 'rng', 42);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

### 2) Advanced simulation with atrophy (2% in left middle frontal gyrus and 3% in right middle frontal gyrus based on Hammers atlas), custom RF field and thicker slices
```matlab
simu = struct('name', 'custom_t1.nii', 'pn', 3, ...
              'resolution', [0.5, 0.5, 1.5], 'rng', []);
simu.atrophy = {'hammers', [28, 29], [2, 3]};
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 3) Thickness simulation (region-wise values, original resolution)
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3, ...
              'resolution', NaN, 'atrophy', [], 'rng', [], ...
              'thickness', [1.5 2.0 2.5]);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

### 4) WMH simulation (medium strength) with simulated RF field
```matlab
simu = struct('name', 'custom_t1.nii', 'pn', 3, 'resolution', NaN, ...
              'WMH', 2, 'rng', []);
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 5) Interactive mode for example 4:
```matlab
simu = struct('pn', 3, 'resolution', NaN, ...
              'WMH', 2, 'rng', []);
rf = struct('percent', 15, 'type', [3, 42], 'save', 0);
mri_simulate(simu, rf);
```

### 6) Apply contrast change (power-law)
```matlab
simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3, ...
              'contrast', 1.3, 'resolution', NaN, 'rng', 0);
rf = struct('percent', 20, 'type', 'A', 'save', 0);
mri_simulate(simu, rf);
```

## File Naming Examples

```
pn3_1.0mm_input_rf20_A.nii                     % Gaussian noise at 3%
snr30_1.0mm_input_rf20_A.nii                   % Rician noise at SNR=30 in WM
pn3_res-10mm_input_thickness2.5mm.nii          % With thickness tag
```
