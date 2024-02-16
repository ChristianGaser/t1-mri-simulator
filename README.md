# mri_simulate
Simulates MR images with various optional features

## Description
This function simulates T1-weighted (T1w) or T2-weighted (T2w) MRI images
based on a high-quality (i.e. 0.5mm spatial resolution) T1w image. By default
Colin27 is used for simulation. The function allows for the inclusion of various 
features in the simulated image such as vessels (only for Colin27), white matter 
hyperintensities (WMHs), Gaussian noise, radiofrequency (RF) B1 inhomogeneities 
(bias field), and options for atrophy or thickness simulation.
If you intend to use any other image for simulation you have to run these steps first:
 * call SPM12 segmentation
 * call CAT12 segmentation in expert mode: save GM, WM, and CSF in native space and deformation 
   field Image->Template (forward)

## Inputs
simu: A structure containing simulation parameters.
 * name: Name of T1w input image.
 * pn: Percentage noise level.
 * rng: Controls the random number generator for consistent noise; empty for default behavior.
 * resolution: Spatial resolution of the simulated image (either a scalar or an xyz-definition
     for anisotropic voxel size).
 * vessel: Boolean flag to add vessels (for T1 image only).
 * WMH: Boolean flag to add white matter hyperintensities.
 * T2: Boolean flag to simulate a T2-weighted image (only possible for Colin27).
 * atrophy: Cell structure that specifies atlas, ROIs and values for simulating atrophy 
     within these ROIs. Multiple ROIs and the respective atrophy values can be defined.
     An atrophy value of 1.5 leads to a GM reduction of about 5, while a value of 2 corresponds to around 10 
     and 3 to around 15. Plese note, that this option takes a lot of time because the atlas labels have to
     be interpolated using categorical interpolation (i.e. each label seperately). Either thickness or atrophy 
     can be simulated.
 * thickness: The WM label of the image is used to add a layer of GM with a defined cortical thickness to have 
     constant thickness. Unlike all other simulations, we can only use the modified label image 
     for the MRI simulation and not the tissue probabilities (which give a more detailed and realistic T1w image).
     Either a scalar value for global constant thickness or a vector with 3 thickness values can be defined
     for the occipital and frontal lobes (1st and 3rd values) and the rest of the brain (2nd value). The 
     Hammer atlas is used to define these areas, as well as subcortical areas and the cerebellum, which are 
     excluded from the thickness simulation to obtain a more realistic MRI. Either thickness or atrophy can be simulated.
 * save: If set to 1 the grounf truth label will be saved.
rf: A structure containing RF bias field parameters.
 * percent: The amplitude of bias field, in percent. Negative values invert the field.
 * type: Specifies the bias field type, options are 'A', 'B', or 'C' for predefined fields from MNI or an
     integer array with two numbers. The first integer value adjusts the local strength of the field by varying 
     maximum frequencies of the FFT. Meaningful values are 1..4, while a value of 3 or 4 corresponds to a bias field
     of a 7T scanner (without further correction such as in mp2rage). The second integer sets the random generator to 
     that seed value, which allows simulating different bias fields.
 * save: If type is a numeric array the simulated bias field can be optionally saved if this value is set to 1.
       

## Deafults
If 'simu' or 'rf' are not provided, they are set to default values:
 simu: {,'name','colin27_t1_tal_hires.nii','pn': 1, 'resolution': 1, 'vessel': 0, 'WMH': 0, 'T2': 0, 'atrophy': {'hammers', [28 29], [1.5 3]}, 'rng', 0}
 rf: {'percent': 20, 'type': [2 0]}

## Output
The function saves the following images:
* simulated MRI image file with the specified features and parameters
* ground truth label with 3 classes
* ground truth PVE label with 3 classes and 2 mixed classes
* bias field

## Usage
mri_simulate(simu, rf);

## Examples
simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3, 'resolution', 0.5, 'vessel', 1, ...
             'WMH', 0, 'T2', 0, 'atrophy', {'hammers', [28 29], [1.5 3]}, 'rng', 0);
rf = struct('percent', 20, 'type', 'A');
mri_simulate(simu, rf);
