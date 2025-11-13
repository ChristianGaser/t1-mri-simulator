# T1-MRI Simulator - Quick Overview

## What is it?

`mri_simulate` is a MATLAB tool that generates realistic simulated T1-weighted MRI brain images with controllable artifacts and anatomical variations. It's designed for testing and validating image analysis algorithms.

---

## Basic Workflow

```
Input T1w Image
      ↓
SPM12 Segmentation → Extract tissue maps (GM/WM/CSF)
      ↓
Apply Modifications (optional):
  • Atrophy (reduce GM in specific regions)
  • Thickness (create uniform cortical thickness)
  • WMH (white matter hyperintensities)
      ↓
Synthesize new T1w image from modified tissues
      ↓
Add RF bias field (optional)
      ↓
Resample to target resolution
      ↓
Add noise (Gaussian or Rician at target WM SNR)
      ↓
Apply contrast change (power-law Y^x, optional)
      ↓
Output: Simulated image + ground truth labels
```

---

## Quick Start

### Minimal Example
```matlab
% Basic simulation with 3% noise
simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

### Common Use Cases

**1. Clean reference image (no artifacts)**
```matlab
simu = struct('name', 'input.nii', 'pn', 0, 'resolution', 1.0);
rf = struct('percent', 0);
mri_simulate(simu, rf);
```

**2. Realistic clinical scan (1.5T typical)**
```matlab
simu = struct('name', 'input.nii', 'pn', 3, 'resolution', 1.0);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

**3. Aging study with atrophy**
```matlab
simu = struct('name', 'input.nii', 'pn', 3, 'resolution', 1.0);
simu.atrophy = {'hammers', [28, 29], [2, 3]};  % 2% and 3% GM loss in ROIs
rf = struct('percent', 20, 'type', 'B');
mri_simulate(simu, rf);
```

**4. White matter disease**
```matlab
simu = struct('name', 'input.nii', 'pn', 3, 'WMH', 2);  % Medium WMH
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```

**5. Cortical thickness study**
```matlab
simu = struct('name', 'input.nii', 'pn', 3, 'thickness', 2.5);  % 2.5mm uniform
rf = struct('percent', 15, 'type', 'A');
mri_simulate(simu, rf);
```

---

## Key Parameters

### Simulation Options (simu)

| Parameter | What it does | Typical values |
|-----------|--------------|----------------|
| `name` | Input T1w image file | Any NIfTI file |
| `pn` | Noise level (% of WM signal) | 0-9 (3 = typical) |
| `resolution` | Output voxel size in mm | NaN, 0.5, 1.0, [1 1 3] |
| `contrast` | Contrast-change | 0.5, 1.5 |
| `WMH` | White matter lesion strength | 0 (off), 1-3 |
| `atrophy` | Regional GM reduction | `{'atlas', [ROIs], [factors]}` |
| `thickness` | Cortical thickness in mm | 1.5-2.5 or `[occ, mid, front]` |
| `rng` | Random seed (reproducibility) | 0 (fixed) or '[]' (random) |

### RF Bias Field Options (rf)

| Parameter | What it does | Typical values |
|-----------|--------------|----------------|
| `percent` | Field strength (%) | ±20 to ±100 |
| `type` | Field pattern | `'A'/'B'/'C'` or `[strength, seed]` |
| `save` | Save bias field file | 0 (no), 1 (yes) |

**RF Field Types:**
- `'A'`, `'B'`, `'C'`: Real MNI-space bias patterns
- `[2, 0]`: Simulated smooth field (strength 2, seed 0)
- `[4, 0]`: Complex 7T-like field (higher strength)

---

## What You Get

Each simulation creates **3-4 output files**:

1. **Simulated image**: `pn3_1.0mm_input_rf20_A.nii`
   - Full brain with all requested effects

2. **Masked image**: `pn3_1.0mm_minput_rf20_A.nii`
   - Brain-only (skull stripped)

3. **Ground truth labels**: `label_pve_1.0mm_input.nii`
   - CSF=1, GM=2, WM=3 (±WMH=4)
   - Useful for training/validation

4. **RF bias field** (if requested): `rf20_A_1.0mm_input.nii`
   - Only for simulated fields
 
 5. **JSON sidecars**: `{simuFile}.json`
       - Metadata with tool info and SimulationParameters (voxel size, pn or snrWM, RF settings, thickness)

---

## Main Features

### 🔬 Tissue Modifications

**Atrophy**
- Reduces gray matter in specific brain regions
- Based on anatomical atlases (Hammers, etc.)
- Increases CSF to simulate volume loss
- Factor 1.5 ≈ 5% reduction, 2.0 ≈ 10%, 3.0 ≈ 15%

**Cortical Thickness**
- Creates uniform cortical ribbon from white matter
- Can vary by region (occipital, middle, frontal)
- Uses distance-based growth with PVE smoothing
- Excludes subcortical structures automatically

**White Matter Hyperintensities (WMH)**
- Simulates age/disease-related white matter changes
- Patchy distribution using random fields
- Based on real WMH probability maps
- Strength 1=mild, 2=moderate, 3=severe

### 🎚️ Image Artifacts

**RF Bias Field**
- Smooth intensity inhomogeneity (B1 field)
- Predefined realistic patterns or custom simulated
- Affects image uniformity across space
- Common in clinical scanners

**Noise & Contrast**
- Gaussian: percentage of white matter mean intensity
- Rician: target SNR in WM (`snrWM`), magnitude noise from complex Gaussian
- Contrast change: power-law mapping Y^x after normalizing to [0,1], then rescaled
- Reproducible with fixed RNG seed

**Resolution Control**
- Isotropic or anisotropic voxels
- Simulates different scanner protocols
- High-quality sinc interpolation

---

## How It Works (Simplified)

1. **Segmentation**: Uses SPM12 to identify GM/WM/CSF in input image

2. **Modification**: Alters tissue distributions based on your parameters
   - Atrophy: reduces GM, increases CSF in ROIs
   - Thickness: grows GM from WM to fixed distance
   - WMH: adds hyperintense patches in white matter

3. **Synthesis**: Recreates T1w image from modified tissue maps
   - Uses original intensity model (Gaussian mixture)
   - Maintains realistic contrast and texture

4. **Artifacts**: Applies bias field and noise to match real scans

5. **Output**: Saves simulated image + ground truth

---

## Requirements

- **MATLAB** (R2017b or later recommended)
- **SPM12** (in MATLAB path)
- **CAT12** (in MATLAB path)
- **Input**: High-quality T1w image (e.g., 0.5mm Colin27 template provided)

---

## Tips & Tricks

### Getting Started
- Use provided `colin27_t1_tal_hires.nii` as template
- Start with default parameters, then customize
- Empty `simu.name` opens file browser (interactive mode)

### Reproducibility
- Set `simu.rng = 0` for identical noise each run
- Set `rf.type = [2, 0]` for fixed bias field
- Document all parameters for your experiments

### Performance
- First run: slow (needs SPM segmentation ~5-10 min)
- Later runs: fast (reuses `*_seg8.mat`)
- Atlas warping adds time (atrophy/thickness options)

### Best Practices
- **Test segmentation first**: Run on input without modifications
- **Validate outputs**: Check if effects match expectations
- **Ground truth**: Use label files for quantitative validation
- **Parameter sweep**: Try multiple noise/bias levels

---

## Common Parameter Combinations

### Study Type Matrix

| Study | Noise (pn) | Resolution | WMH | Atrophy | Thickness | RF Field |
|-------|------------|------------|-----|---------|-----------|----------|
| **Algorithm test** | 0 | Original | 0 | No | No | 0% |
| **Clinical 1.5T** | 3 | 1.0mm | 0-1 | No | No | 20%, A/B/C |
| **Clinical 3T** | 2 | 0.8mm | 0-1 | No | No | 20%, A/B/C |
| **7T research** | 2 | 0.7mm | 0 | No | No | 40%, [4,0] |
| **Aging study** | 3 | 1.0mm | 1-3 | Yes | No | 20%, A/B/C |

---

## Limitations

- **Single modality**: Only T1w (no T2, FLAIR, etc.)
- **No motion**: Motion artifacts not yet implemented
- **Simplified WMH**: Single intensity class
- **Requires segmentation**: Input must be segmentable by SPM

---

## Output Interpretation

### File Naming
```
pn{noise}_{resolution}mm_{input}{options}.nii  or  snr{SNR}_{resolution}mm_{input}{options}.nii

Examples:
pn3_1.0mm_colin27_t1_tal_hires.nii                    → Basic
pn3_0.5mm_input_rf20_A.nii                            → With RF field A
pn3_1.0mm_input_rf15_2_42_WMH2.nii                    → Simulated RF + WMH
pn3_1.0mm_input_hammers_28_2.nii                      → Atrophy in ROI 28
pn3_1.0mm_input_thickness2.5mm.nii                    → Uniform 2.5mm cortex
pn3_1.0mm_input_thickness1.5mm-2.5mm.nii             → Regional thickness
snr30_1.0mm_input_rf20_A.nii                          → Rician noise at WM SNR=30
```

### Label Values
- **1.0** = Pure CSF
- **2.0** = Pure GM
- **3.0** = Pure WM
- **4.0** = WMH (if enabled)
- **Intermediate** = Partial volume (e.g., 2.5 = 50% GM + 50% WM)

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "seg8.mat not found" | Normal on first run - SPM will create it |
| "Atlas not found" | Check atlas name spelling |
| "Cannot combine atrophy and thickness" | Choose only one |
| Simulation too slow | Atlas warping takes time; use smaller ROIs or skip |
| Unrealistic results | Check input image quality and parameters |

---

## Interactive Mode

Don't know which file to use? Just run:
```matlab
simu = struct('pn', 3);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
```
A file browser will open, and you can select any T1w image!

---

## Further Reading

- **Detailed documentation**: See `T1_SIMULATION_SCHEME.md` for technical details
- **Code comments**: `mri_simulate.m` has extensive inline documentation
- **SPM12 manual**: https://www.fil.ion.ucl.ac.uk/spm/doc/
- **CAT12 manual**: http://www.neuro.uni-jena.de/cat12-help/

---

## Quick Reference Card

### Minimal Working Examples

```matlab
% Default everything
mri_simulate();

% Just noise
simu = struct('name', 'input.nii', 'pn', 5);
mri_simulate(simu);

% Just bias field
simu = struct('name', 'input.nii', 'pn', 0);
rf = struct('percent', 30, 'type', 'B');
mri_simulate(simu, rf);

% Everything combined
simu = struct('name', 'input.nii', 'pn', 3, 'resolution', 1.0, 'WMH', 2);
simu.atrophy = {'hammers', [28], [2]};
rf = struct('percent', 20, 'type', [3, 42], 'save', 1);
mri_simulate(simu, rf);
```

---

**Version**: 1.0  
**Author**: Christian Gaser  
**Repository**: github.com/ChristianGaser/t1-mri-simulator  
**License**: See LICENSE file
