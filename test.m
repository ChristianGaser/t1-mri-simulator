function mri_simulate(simu, rf)
% MRI_SIMULATE - Simulates magnetic resonance imaging (MRI) with customizable features.
%
% Overview:
%   `mri_simulate` generates simulated MRI images, allowing for T1-weighted (T1w) or T2-weighted (T2w) imaging
%   simulations. It employs a high-resolution (e.g., 0.5mm) T1w image, typically Colin27, as a base. Users can introduce
%   various artifacts and features like vessels, white matter hyperintensities (WMHs), Gaussian noise, and RF B1 inhomogeneities.
%   It supports simulations of atrophy or cortical thickness modifications. Preprocessing with SPM12 and CAT12 segmentation
%   is required for custom images.
%
% Syntax:
%   mri_simulate(simu, rf)
%
% Parameters:
%   simu (struct): Simulation parameters.
%       - 'name' (string): Filename of the T1w input image.
%       - 'pn' (double): Percentage noise level to introduce Gaussian noise.
%       - 'rng' (double or []): Seed for the random number generator; use [] for MATLAB's default behavior.
%       - 'resolution' (double or [x, y, z]): Spatial resolution of the simulated image.
%       - 'vessel' (logical): Flag to add vessels to the T1 image. Only applicable for Colin27.
%       - 'WMH' (logical): Flag to add white matter hyperintensities.
%       - 'T2' (logical): Flag for simulating a T2-weighted image. Only works with Colin27.
%       - 'atrophy' (cell): Specifies regions of interest (ROIs) for simulating atrophy, including atlas name, ROI IDs, and atrophy values.
%       - 'thickness' (double or [double]): Specifies the cortical thickness for simulation.
%       - 'save' (logical): Flag to save the ground truth label if set to 1.
%   rf (struct): RF bias field parameters.
%       - 'percent' (double): Amplitude of the bias field in percentage. Negative values invert the field.
%       - 'type' (char or [int, int]): Specifies the bias field type or an array defining field strength and seed for random variation.
%       - 'save' (logical): Option to save the simulated bias field if 'type' is numeric and 'save' is set to 1.
%
% Optional Inputs:
%   Default values are used if 'simu' or 'rf' parameters are not provided. See examples for default structures.
%
% Outputs:
%   Simulated MRI image file based on the specified parameters and features.
%
% Usage:
%   To simulate an MRI, specify the simulation (`simu`) and RF bias field (`rf`) parameters:
%       mri_simulate(simu, rf);
%
% Examples:
%   Example 1 - Basic simulation with added vessels and specific noise:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3, 'resolution', 0.5, 'vessel', true, ...
%                     'WMH', false, 'T2', false, 'atrophy', {}, 'rng', 42);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%   Example 2 - Advanced simulation with atrophy and custom RF field:
%       simu = struct('name', 'custom_t1.nii', 'pn', 2, 'resolution', [0.5, 0.5, 0.5], 'vessel', false, ...
%                     'WMH', true, 'T2', false, 'atrophy', {'hammers', [28, 29], [2, 3]}, 'rng', []);
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
% Notes:
%   - Preprocessing steps are required for using custom images for simulation. Refer to SPM12 and CAT12 documentation.
%   - Simulating atrophy or cortical thickness requires significant computational resources due to the interpolation of atlas labels.
%   - Ensure the MATLAB path includes dependencies such as SPM12 and CAT12 toolboxes for processing custom images.
%
% See also: SPM12, CAT12
