function mri_simulate(simu, rf)
% MRI_SIMULATE - Simulates MR images with various optional features
%
% Overview:
%   `mri_simulate` generates simulated MRI images, allowing for T1-weighted
%   (T1w) imaging simulations. It employs a high-quality (high-resolution, low
%   noise) T1w image as a base. Users can introduce various artifacts and 
%   features like white matter hyperintensities (WMHs), Gaussian noise, and 
%   RF B1 inhomogeneities. 
%   It supports simulations of atrophy or cortical thickness modifications. 
%   Preprocessing with SPM12 segmentation is required for custom images.
%
%   Thickness/PVE pipeline (when simu.thickness is set):
%   - A constant cortical thickness (global or region-wise) is synthesized by
%     expanding GM outward from the original WM using a Euclidean distance map.
%   - To obtain partial-volume-like transitions, the label boundary is shifted
%     across 20 sub-voxel offsets in the range [-0.25, 0.25] voxels, simulating
%     realistic boundary uncertainty. Each offset yields a hard label image
%     (CSF=1, GM=2, WM=3), and the results are averaged to produce a smooth
%     PVE-like label map.
%   - The final PVE label map is then used to synthesize a T1 image from the
%     SPM segmentation model by replacing GM/WM/CSF class posteriors and using
%     their Gaussian mixture parameters (means, variances, weights).
%
% Syntax:
%   mri_simulate(simu, rf)
%
% Parameters:
%   simu (struct): Simulation parameters. Defaults are applied for missing fields.
%       - 'name' (char): T1-weighted input image filename. Default: '' (empty),
%         which triggers an interactive file selection dialog (T1 only).
%       - 'pn' (double): Percentage noise level (Gaussian) relative to the WM
%         mean intensity. Ignored if 'snrWM' > 0. Default: 3 (percent of WM).
%       - 'snrWM' (double): If >0, adds Rician noise at a user-defined SNR for
%         white matter. Uses the (noise-free) WM mean to compute the complex
%         noise sigma via sigma = WMmean / snrWM, and generates magnitude
%         Rician noise: sqrt((S + n1).^2 + n2.^2). Default: 0 (disabled).
%       - 'rng' (double or []): Seed for the random number generator. Default: 0
%         (reproducible noise across runs). Set [] to use MATLAB's default RNG behavior
%         (non-deterministic across sessions).
%       - 'resolution' (double or [x, y, z]): Spatial resolution of the
%         simulated image. Default: NaN (keep original resolution). If scalar,
%         it is applied to all three axes; if a 3-vector, each axis is set individually.
%       - 'WMH' (integer or >=1 scalar): Strength of simulated white matter
%         hyperintensities (WMHs).
%           0  -> no WMHs
%           1  -> mild, mainly periventricular patches
%           2  -> medium extent/contrast
%           3  -> strong, widespread WMHs
%         Notes:
%           • Values >=1 are allowed (not just 1/2/3). Higher values broaden
%             the WMH extent by reshaping the WMH prior with an exponent
%             1/(WMH-0.8), and adjust the contribution to the label map by a
%             scaling of ~1/WMH^0.75 to keep intensities in a plausible range.
%           • WMHs are constrained to (eroded) WM and modulated by a random
%             field to produce patchy distributions.
%         Default: 0 (disabled).
%       - 'atrophy' (cell): Specifies regions of interest (ROIs) for simulating
%         atrophy, including atlas name, ROI IDs, and atrophy values. Multiple ROIs 
%         and the respective atrophy values can be defined. An atrophy value of 
%         1.5 leads to a GM reduction of about 5%, while a value of 2 corresponds 
%         to around 10% and 3 to around 15%. Please note, that this option takes 
%         a lot of time because the atlas labels have to be interpolated using 
%         categorical interpolation (i.e. each label separetely). Either thickness 
%         or atrophy can be simulated. Default: [] (disabled).
%       - 'thickness' (double or [double double double]): Specifies the cortical 
%         thickness for simulation. The WM label of the image is used to add a 
%         layer of GM with a defined cortical thickness to have constant thickness. 
%         Unlike all other simulations, we can only use the modified label image 
%         for the MRI simulation and not the tissue probabilities (which give a 
%         more detailed and realistic T1w image).
%         Either a scalar value for global constant thickness or a vector with 
%         3 thickness values can be defined for the occipital and frontal lobes 
%         (1st and 3rd values) and the rest of the brain (2nd value). The 
%         Hammer atlas is used to define these areas, as well as subcortical areas 
%         and the cerebellum, which are excluded from the thickness simulation to 
%         obtain a more realistic MRI. Either thickness or atrophy can be simulated.
%         Default: 0 (disabled).
%   rf (struct): RF bias field parameters.
%       - 'percent' (double): Amplitude of the bias field in percentage.
%         Negative values invert the field. Default: 20.
%       - 'type' (char or [int, int]): Specifies the bias field type:
%           'A' | 'B' | 'C' for predefined MNI fields, or [strength, rngSeed]
%           for a simulated field. Meaningful strength values are 1..4; 3–4
%           resemble stronger inhomogeneities (e.g., 7T without correction).
%         Default: [2 0] (simulated field with moderate strength and RNG seed 0).
%       - 'save' (logical): Save the simulated bias field only when 'type' is
%         numeric (simulated). Ignored for predefined 'A'/'B'/'C'. Default: 0.
%
% Optional Inputs and Defaults:
%   - If 'simu' and/or 'rf' are omitted, default values are applied internally.
%     You can also pass a partial struct; any missing fields are filled with
%     defaults (i.e., user-specified fields override defaults only where set).
%   - Interactive mode: if 'simu.name' is empty (''), the function opens a file
%     selection dialog (spm_select) and runs the simulation for the chosen file(s).
%     This allows you to launch the tool without pre-specifying the input path.
%   See the examples below for constructing minimal 'simu' and 'rf' structs.
%
% Outputs:
%   Simulated MRI image file based on the specified parameters and features.
%
% Usage:
%   To simulate an MRI, specify the simulation (`simu`) and RF bias field
%   (`rf`) parameters:
%       mri_simulate(simu, rf);
%
% Examples:
%   Example 1 - Basic simulation with specific noise with 0.5mm voxel size:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3,...
%                     'resolution', 0.5,...
%                     'atrophy', [], 'rng', 42);
%       rf = struct('percent', 20, 'type', 'A','save',0);
%       mri_simulate(simu, rf);
%
%   Example 2 - Advanced simulation with atrophy (2% in left middle frontal gyrus 
%               and 3% in right middle frontal gyrus based on Hammers atlas), 
%               custom RF field and thicker slices: 
%       simu = struct('name', 'custom_t1.nii', 'pn', 3,...
%                     'resolution', [0.5, 0.5, 1.5],...
%                     'rng', []);
%       simu.atrophy = {'hammers',[28, 29], [2, 3]};
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%   Example 3 - Thickness simulation with 3 different thickness values using 
%   original voxel size:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3,...
%                     'resolution', NaN,...
%                     'atrophy', [], 'rng', [],...
%                     'thickness', [1.5 2.0 2.5]);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
%   Example 4 - Simulation with custom RF field and added WMHs (medium strength)
%       simu = struct('name', 'custom_t1.nii', 'pn', 3,...
%                     'resolution', NaN, 'WMH', 2, ...
%                     'rng', []);
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%
%
% TODO: simulation of motion artefacts using FFT and shift of phase information

% Default simulation parameters
def.name       = '';
def.pn         = 3;
def.resolution = NaN;
def.WMH        = 0;
def.atrophy    = [];
def.thickness  = 0;
def.rng        = 0;
def.snrWM      = 0;

if nargin < 1, simu = def;
else, simu = cat_io_checkinopt(simu, def); end

% Default bias field parameters
def.percent    = 20;
def.type       = [2 0];
def.save       = 0;

if nargin < 2, rf = def;
else, rf = cat_io_checkinopt(rf, def); end

% call tool interactively
if isempty(simu.name)
  simu.name = spm_select(Inf,'image','Select T1w image(s) for simulation');
end

n_images = size(simu.name,1);
if n_images > 1
  P = simu.name;
  for i=1:n_images
    simu.name = deblank(P(i,:));
    mri_simulate(simu, rf);
  end
  return
end

% iterations for correction for regions with too large thickness values
% not working properly!!!
n_thickness_corrections = 0;

[pth, name, ext, ~] = spm_fileparts(simu.name);
if strcmp(ext,'.gz')
  fname = gunzip(fullfile(pth, [name ext]));
  simu.name = fname{1};
  [pth, name, ext] = spm_fileparts(simu.name);
  is_gz = 1;
else
  is_gz = 0;
end
  
pth_root = fileparts(which(mfilename));

% name of seg8.mat file that contains SPM12 segmentation parameters
mat_name = fullfile(pth, [name '_seg8.mat']);

% CAT12 template dir is later used for defining atrophy atlas
template_dir = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym');

% call SPM segmentation if necessary and only save the seg8.mat file
if ~exist(mat_name,'file')
  matlabbatch{1}.spm.spatial.preproc.channel.vols = {simu.name};
  spm_jobman('run',matlabbatch);
  clear matlabbatch

  % remove native segmentations that were saved
  for i=1:5
    spm_unlink(fullfile(pth,['c' num2str(i) name ext]));
  end
end

% check that only one of these options is used
if ((isfield(simu,'atrophy') && ~isempty(simu.atrophy) && any(simu.atrophy{3} > 1))) && any(simu.thickness)
  fprintf('Option to simulate atrophy cannot be combined with option to simulate thickness images.\n');
  return
end

% check that strength for simulated bias field is at least 1
if isnumeric(rf.type) && rf.type(1) < 1
  fprintf('Strength of simulated bias field should be at least 1.\n');
  return
end

% check that rf.save is not set for predefined MNI bias fields
if ischar(rf.type) && rf.save
  fprintf('Predifined MNI bias fields cannot be saved.\n');
  rf.save = 0;
end

% any parameter to simulate atrophy defined?
if isfield(simu,'atrophy') && ~isempty(simu.atrophy) && any(simu.atrophy{3} > 1)
  simu_atrophy = 1;
  atlas_name = fullfile(template_dir,[simu.atrophy{1} '.nii']);
  
  % check that atlas exists
  if ~exist(atlas_name,'file')
    fprintf('Atlas %s could not be found. Please check that the name is correct.\n', atlas_name);
    return
  end
else
  simu_atrophy = 0;
end

if simu.WMH < 1 && simu.WMH ~= 0
  fprintf('Strength of simulating WMHs should be either 0 or >=1.\n');
  return
end

% load seg8.mat file and define some parameters
res = load(mat_name);

if size(res.mn,1) > 1
  fprintf('Multi-modal segmentation is not recommended! Try again to segment the image using T1w data only.\n');
end

% get means for GM/WM/CSF
mn = zeros(3,1);
for k=1:3
  ind = find(res.lkp==k);
  for j=1:numel(ind)
    mn(k) = mn(k) + res.mg(ind(j))*res.mn(1,ind(j));
  end  
end

% check that it's indeed T1w data by checking CSF < GM < WM
[~, ind] = sort(mn);
if ind ~= [3 1 2]
  fprintf('Warning: No typical T1w intensities were found. Please note that segmentation quality can be much lower for non T1w data.\n');
  if simu.WMH
    simu.WMH = 0;
    fprintf('Warning: WMHs can be only simulated for T1w data. WMH option was therefore disabled.\n');
  end
end

[~, bname, ext] = spm_fileparts(res.image(1).fname);
res.image(1) = spm_vol(fullfile(pth,[bname ext]));
idef_name = fullfile(pth, ['iy_' name ext]);
V = res.image(1);
dim   = V.dim(1:3);
vx = sqrt(sum(V.mat(1:3,1:3).^2));

% obtain SPM segmentations and write inverse deformation field
[Ysrc, Ycls] = cat_spm_preproc_write8(res,zeros(max(res.lkp),4),zeros(2,2),[1 0],0,2);

% get tissue thresholds for CSF/GM/WM
T3th = get_tissue_thresholds(Ysrc, Ycls, res);

% LAS correction and SANLM denoising
Ycorr = cat_main_LASsimple(Ysrc, Ycls, T3th);
cat_sanlm(Ycorr,3,1);

% Replace GM/WM/CSF segmentation by labels using LAS corrected image
Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
Yseg = zeros([dim, 3]);

% Use CAT12 adaptive probability region-growing (APRG) approach for
% skull-stripping
brainmask = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th);
Ycorr = Ycorr.*brainmask;

seg_order = [2 3 1];
for i = 1:3
    Yseg(:,:,:,i) = Yp0toC(3*Ycorr, seg_order(i)); 
end
clear Ycls

% add atrophy to GM by decreasing GM value in ROI and increasing CSF value
if simu_atrophy
  Yseg = simulate_atrophy(simu, Yseg, dim, template_dir, idef_name);
end

% get ground truth label using GM/WM/CSF
label_pve = zeros(dim, 'single');
order = [3 1 2];
for k = 1:3
    label_pve = label_pve + k*(Yseg(:,:,:,order(k)));
end

% load WMH map
if simu.WMH
  [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name);
else
  WMH = [];
end

% extend target voxel size if needed
if isscalar(simu.resolution)
  if isnan(simu.resolution)
    simu.resolution = vx;
    change_resolution = 0;
  else
    simu.resolution = simu.resolution * ones(1,3);
    change_resolution = 1;
  end
else
  if isnan(simu.resolution(1))
    simu.resolution = vx;
    change_resolution = 0;
  else
    change_resolution = 1;
  end
end

if any(simu.thickness)
  [label_pve, Yseg] = simulate_thickness(label_pve, simu, Yseg, dim, ...
        template_dir, idef_name, vx, order, n_thickness_corrections);
end

Ysimu = synthesize_from_segmentation(Yseg, name, res, mn, dim, WMH);

% apply either predefined MNI bias field or simulated bias filed before resizing
% to defined output resolution
if rf.percent ~= 0
  if ischar(rf.type)
    [Ysimu, rf_field] = add_bias_field(Ysimu, rf, idef_name, pth_root); % add predefined MNI field
  else
    [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf);
  end
end

mx_vol = max(Ysimu(:));

% output matrix
Vres.dim = round(V.dim.*vx./simu.resolution);
P = spm_imatrix(V.mat);
P(7:9) = P(7:9)./vx.*simu.resolution;
P(1:3) = P(1:3) + vx - simu.resolution;
Vres.mat = spm_matrix(P);

% output in defined resolution
volres   = zeros(Vres.dim);
labelres_pve = zeros(Vres.dim);

if change_resolution
  for sl = 1:Vres.dim(3)
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    M1 = Vres.mat\V.mat\M;
    
    % use sinc interpolation for simulated image
    volres(:,:,sl) = spm_slice_vol(Ysimu,M1,Vres.dim(1:2),-5);
    % and linear interpolation for label image
    labelres_pve(:,:,sl) = spm_slice_vol(label_pve,M1,Vres.dim(1:2),1);
  end
else % we can skip interpolation if voxels size is the same
  volres = Ysimu;
  labelres_pve = label_pve;
end

volres = volres / mx_vol;

% optionally ensure same noise for every trial
if ~isempty(simu.rng)
  rng(simu.rng,'twister');
end

% Add noise:
% - If simu.snrWM>0: add Rician noise with target WM SNR
% - Else: fall back to Gaussian noise using percentage of WM mean
if isfield(simu,'snrWM') && ~isempty(simu.snrWM) && simu.snrWM > 0
  % Desired SNR is defined for the underlying complex signal amplitude at WM
  % Compute complex noise std using absolute WM mean, then convert to normalized units
  sigma_abs = mn(2) / simu.snrWM;    % absolute-domain sigma (same units as Ysimu)
  sigma = sigma_abs / mx_vol;        % normalized-domain sigma
  n1 = sigma * randn(size(volres));
  n2 = sigma * randn(size(volres));
  volres = sqrt( (volres + n1).^2 + (n2).^2 );
  % Clamp and rescale back
  volres = mx_vol * max(min(volres, 1), 0);
else
  % Gaussian noise using percentage of WM mean intensity
  sigma_abs = (simu.pn/100) * mn(2); % absolute-domain std dev
  sigma = sigma_abs / mx_vol;        % normalized-domain std dev
  noise = sigma * randn(size(volres));
  volres = volres + noise;
  volres = mx_vol * max(min(volres, 1), 0);
end

str = '';
if (exist('rf','var') && rf.percent ~= 0)
  if isnumeric(rf.type)
    str = sprintf('_rf%g_%d_%d', rf.percent, rf.type(1), rf.type(2));
  else
    str = sprintf('_rf%g_%s', rf.percent, rf.type);
  end
end
if simu.WMH
  str = sprintf('%s_WMH%g', str, simu.WMH);
end
if simu_atrophy
  if numel(simu.atrophy{2}) > 1
    str = sprintf('%s_%s_multi',str,simu.atrophy{1});
  else
    str = sprintf('%s_%s_%d_%g',str,simu.atrophy{1},simu.atrophy{2},simu.atrophy{3});
  end
end
if any(simu.thickness)
  if isscalar(simu.thickness)
    str = sprintf('%s_thickness%gmm',str,simu.thickness);
  else
    str = sprintf('%s_thickness%gmm-%gmm',str,min(simu.thickness),max(simu.thickness));
  end
end

mean_resolution = round(10*mean(simu.resolution))/10;

% write simulated image
simu_name = fullfile(pth, sprintf('pn%g_%gmm_%s%s.nii',simu.pn, mean_resolution, name, str));
fprintf('Save simulated image %s\n', simu_name);
Vres.fname = simu_name;
Vres.pinfo = [1 0 352]';
Vres.dt    = [16 0];
spm_write_vol(Vres, volres);

% write simulated masked image
simu_name = fullfile(pth, sprintf('pn%g_%gmm_m%s%s.nii',simu.pn, mean_resolution, name, str));
fprintf('Save simulated masked image %s\n', simu_name);
Vres.fname = simu_name;
mind = labelres_pve(:,:,:) > 0.5;
spm_write_vol(Vres, volres.*mind);

% write ground truth label
if simu.WMH
  WMHstr = sprintf('_WMH%g', simu.WMH);
else
  WMHstr = '';
end
if any(simu.thickness)
  label_pve_name = fullfile(pth, sprintf('label_pve_%gmm_%s%s%s.nii', mean_resolution, name, str, WMHstr));
else
  label_pve_name = fullfile(pth, sprintf('label_pve_%gmm_%s%s.nii', mean_resolution, name, WMHstr));
end

fprintf('Save %s\n', label_pve_name);
Vres.fname = label_pve_name;
Vres.pinfo = [1/255/3 0 352]';
Vres.dt    = [4 0];
spm_write_vol(Vres, labelres_pve);

% save simulated bias field if defined
if rf.save
  rf_name = fullfile(pth, sprintf('%s_%gmm_%s.nii', str(2:end), mean_resolution, name));
  fprintf('Save %s\n', rf_name);
  Vres.fname = rf_name;
  Vres.pinfo = [1/max(rf_field(:)) 0 352]';
  Vres.dt    = [2 0];
  spm_write_vol(Vres, rf_field);
end

% remove temporary files
spm_unlink(idef_name);
if is_gz
  spm_unlink(simu.name);
end

fprintf('================================================================================\n');

%==========================================================================
% function t = transf(B1,B2,B3,T)
% from spm_preproc_write8.m
%
% Purpose
%   Reconstruct a low-rank 3D field (e.g., bias field) using separable DCT
%   bases along x/y/z with coefficients T.
%
% Inputs
%   B1,B2,B3 - DCT basis matrices for x, y, and z.
%   T        - Coefficient tensor; if empty, a zero field is returned.
%
% Output
%   t        - Reconstructed field with size [size(B1,1) size(B2,1) size(B3,1)].
%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
  d2 = [size(T) 1];
  t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
  t  = B1*t1*B2';
else
  t = zeros(size(B1,1),size(B2,1),size(B3,1));
end
return


%==========================================================================
% function p = likelihoods(f,bf,mg,mn,vr)
% from spm_preproc_write8.m
%
% Purpose
%   Compute per-voxel likelihoods for a Gaussian mixture model given features
%   (optionally bias-corrected) and component parameters.
%
% Inputs
%   f   - cell array with one feature image; values are vectorized internally.
%   bf  - optional bias field multiplier (same shape as f{1}); [] for none.
%   mg  - mixture weights for each Gaussian component (1 x K).
%   mn  - means per component (N x K, N=#features).
%   vr  - covariance matrices per component (N x N x K).
%
% Output
%   p   - per-voxel likelihoods for each component (numel(f{1}) x K).
%==========================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N
  if isempty(bf)
    cr(:,n) = double(f{1}(:));
  else
    cr(:,n) = double(f{1}(:).*bf{1}(:));
  end
end

p  = ones(numel(f{1}),K);
for k=1:K
  amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
  d      = bsxfun(@minus,cr,mn(:,k)')/chol(vr(:,:,k));
  p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
return


%==========================================================================
% function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, order, n_thickness_corrections)
%
% Purpose
%   Synthesize a constant cortical thickness by growing GM outward from the
%   original WM using an Euclidean distance map, and convert the resulting
%   hard labels into a partial-volume-like (PVE) segmentation by boundary
%   jittering and averaging.
%
% Inputs
%   label    - single(dims): Current PVE-like label image with values in [1..3]
%              (CSF=1, GM=2, WM=3). Used as baseline and to preserve subcortical
%              and cerebellar regions.
%   simu     - struct: Simulation options. Relevant fields:
%                .thickness: either a scalar (global thickness, in mm) or a
%                            3-element vector [occipital rest frontal], in mm.
%   Yseg     - single(dims,3): Current tissue probability/label volumes in the
%              order specified by 'order'. Will be overwritten by the simulated
%              PVE-like segmentation created here.
%   d        - [nx ny nz]: Volume dimensions.
%   template_dir - char: Path to CAT/SPM template directory (for Hammer atlas).
%   idef_name    - char: Filename of inverse deformation field to warp atlas
%                        into subject/native space with categorical interpolation.
%   vx       - [vx vy vz]: Voxel sizes in mm.
%   order    - [3x1] int: Mapping from class index (CSF/GM/WM) to Yseg order.
%   n_thickness_corrections - int: Optional iterations for correcting regions
%                        where the requested thickness is too large for sulci;
%                        set to 0 for default behavior.
%
% Outputs
%   label    - single(dims): New PVE-like label map in [1..3] after averaging
%              across boundary jitters, aligned with Yseg/order.
%   Yseg     - single(dims,3): Updated class volumes (CSF/GM/WM) representing
%              the simulated PVE-like segmentation.
%
% Algorithm
%   1) Atlas masks: Warp Hammer atlas to native space and build masks to exclude
%      subcortical/cerebellar regions from thickness manipulation. Optionally
%      define region masks to apply three distinct thickness values (occipital,
%      rest, frontal) when simu.thickness is a 3-vector.
%   2) Boundary jittering (PVE simulation): To emulate partial volume effects,
%      shift the label boundaries by 20 sub-voxel offsets uniformly spaced in
%      [-0.25, 0.25] voxels. For each offset:
%        a. Threshold labels to obtain a hard WM mask and clean it with simple
%           morphological steps.
%        b. Compute Euclidean distance transform from WM (compensated by 0.5*voxel).
%        c. Set CSF everywhere, then assign GM to voxels where D_WM <= thickness
%           (per region), preserving WM where present.
%        d. Optionally correct gyri too thin for the requested thickness by
%           widening sulci (n_thickness_corrections>0).
%        e. Convert the hard labels (1..3) to one-hot class volumes and add to
%           an accumulator.
%   3) Average the 20 accumulated volumes to form a smooth PVE-like segmentation
%      and rebuild the final label image in [1..3] by weighted sum with class IDs.
%
% Notes
%   - Class encoding: CSF=1, GM=2, WM=3 throughout.
%   - When simu.thickness is scalar, the same thickness is applied globally.
%     When it has 3 values, they are applied to occipital, rest, and frontal
%     regions as defined by the Hammer atlas masks.
%   - This function modifies Yseg directly to reflect the new PVE-like tissue
%     maps, which are later used by the synthesis step to generate a T1 image.
%==========================================================================
function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, order, n_thickness_corrections)

Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
csf_val = 1; gm_val = 2; wm_val = 3;

% warp atlas to native space using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
hammers_name = fullfile(template_dir,'hammers.nii');
hammers = cat_vol_defs(struct('field1',{{idef_name}},'images',{{hammers_name}},'interp',-1,'modulate',0));
hammers = single(hammers{1}{1});

% create mask for basal ganglia and cerebellum
mask_orig = (hammers>0 & hammers<5) | hammers==9 | hammers==10 | ...
            (hammers>16 & hammers<20) | (hammers>33 & hammers<50) | ...
             hammers==74 | hammers==75;
mask_orig = cat_vol_morph(mask_orig,'dd',2);

% create mask for occipital and frontal areas and remaining parts
mask_thickness{1} = (hammers > 63.5 & hammers < 67.5) | (hammers > 21.5 & hammers < 23.5); % occipital
mask_thickness{3} = (hammers > 55.5 & hammers < 59.5) | (hammers > 27.5 & hammers < 29.5) | (hammers > 68.5 & hammers < 71.5); % frontal
mask_thickness{2} = ~mask_thickness{1} & ~mask_thickness{3}; % remaining parts

mask = round(label) > 0;

% force stronger PVE effects by smoothing
spm_smooth(label,label,2.5*vx);

label0 = label;

label1 = cell(numel(simu.thickness),1);

Yseg(:,:,:,1:3) = 0;

% vary range of PVE from -0.25..0.25 in 20 steps
pve_range = linspace(-0.25,0.25,20);

for pve_step = 1:numel(pve_range)

  % define wm, remove disconnected regions and dilate it to thickne thin gyri
  wm  = round(label+pve_range(pve_step)) == wm_val;
  wm = cat_vol_morph(wm,'l',1);
  wm = cat_vol_morph(wm,'dc',1);
  
  % euclidean distance to wm
  Dwm = (bwdist(wm) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance

  for k=1:numel(simu.thickness)
  
    label1{k} = round(label+pve_range(pve_step));

    label1{k}(~wm) = csf_val;
    label1{k}(~mask) = 0;
    
    % limit dilated gm to defined thickness
    label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k)) = gm_val;  

    for i=1:n_thickness_corrections
      % check distance inside gm
      D = (bwdist(label1{k}==csf_val) + bwdist(label1{k}==wm_val) - 1)*mean(vx);
      D(label1{k}~=gm_val) = 0;
          
      % dilate gm where distance in gm is too large to make more space
      gm1 = cat_vol_morph(D > 0.97*simu.thickness(k),'o',1);
      gm1 = cat_vol_morph(gm1,'dd',1);
      gm1(~mask) = 0;
      
      % widen sulci is thickness was too large for that area
      label1{k}(gm1 > 0 | label1{k} == gm_val) = csf_val;
      
      % euclidean distance to wm
      wm  = label1{k}==wm_val;
      Dwm = (bwdist(wm) - 0.5) * mean(vx);
      
      % limit dilated gm to defined thickness
      label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k)) = gm_val;  
    end
    
  end
    
  % replace tissue maps with modified label
  for j = 1:3
  
    if isscalar(simu.thickness)
      % only simulate 2mm thickness
      tmp_seg = single(round(label1{1}) == (j));
    else
      tmp_seg = Yseg(:,:,:,order(j));
      for k = 1:3
        tmp_seg(mask_thickness{k}) = single(round(label1{k}(mask_thickness{k})) == (j));
      end
    end
    % get original label for basal ganglia and cerebellum
    tmp_seg(mask_orig) = Yp0toC(label0(mask_orig), j);

    Yseg(:,:,:,order(j)) = Yseg(:,:,:,order(j)) + tmp_seg/numel(pve_range);
  end
end

% update ground truth label
label = zeros(d, 'single');

for k = 1:3
    label = label + k*Yseg(:,:,:,order(k));
end

%==========================================================================
% function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name)
%
% Purpose
%   Apply regional atrophy by increasing CSF (and effectively reducing GM)
%   within specified atlas ROIs warped into native space.
%
% Inputs
%   simu         - struct with .atrophy = {atlasName, roiIds[], factors[]}
%   Yseg         - single(dims,3): tissue maps (CSF/GM/WM order as used).
%   dims         - [nx ny nz] dimensions of the volume.
%   template_dir - path to atlas templates; atlas is warped categorically.
%   idef_name    - inverse deformation field to native space.
%
% Output
%   Yseg         - updated tissue maps with CSF increased in target ROIs.
%==========================================================================
function Yseg = simulate_atrophy(simu, Yseg, dims, template_dir, idef_name)

% warp atlas to native spave using categorical interpolation
fprintf('Transform atlas to native space. This may take a while...\n');
atlas_name = fullfile(template_dir,[simu.atrophy{1} '.nii']);
atlas = cat_vol_defs(struct('field1',{{idef_name}},'images',{{atlas_name}},'interp',-1,'modulate',0));
atlas = single(atlas{1}{1});

% go through all defined ROIs
for i = 1:numel(simu.atrophy{2})
  ind_atlas = round(atlas) == simu.atrophy{2}(i);
  
  % check that ROI is not empty because it was not correct defined
  if isempty(ind_atlas)
    fprintf('ROI #%d seems to be wrong and does not exist in %s.\n', simu.atrophy{2}, simu.atrophy{1});
    return
  end
  
  % create atrophy mask with defined value in mask and otherwise 1
  mod_atlas = zeros(dims, 'single');
  mod_atlas(ind_atlas) = simu.atrophy{3}(i);
  mod_atlas(~ind_atlas) = 1.0;
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas = sum(maxind(ind_atlas) == 1);

  % increase CSF by factor defined in atrophy mask
  Yseg(:,:,:,3) = mod_atlas.*Yseg(:,:,:,3);
  
  [~, maxind] = max(Yseg,[],4);
  GM_atlas_simu = sum(maxind(ind_atlas) == 1);
  fprintf('GM reduction in region #%d = %3.2f%s\n',simu.atrophy{2}(i),100*(1-GM_atlas_simu/GM_atlas),'%');
end


%==========================================================================
% function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d, WMH)
%
% Purpose
%   Generate a T1-like image from provided CSF/GM/WM maps by inserting them
%   into SPM's mixture model and computing the expected intensity per voxel.
%
% Inputs
%   vol_seg - single(dims,3): tissue maps used in place of SPM posteriors.
%   name    - base name for progress display.
%   res     - struct from SPM segmentation (mg, mn, vr, Tbias, etc.).
%   mn      - 3x1 means for CSF/GM/WM replacing the corresponding mixture means.
%   d       - [nx ny nz] dimensions.
%   WMH     - optionally added white matter hyperintensities (WMHs)
%
% Output
%   Ysimu   - synthesized T1-weighted image volume (single).
%==========================================================================
function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d, WMH)
% go through all peaks that are defined
% mainly copied from spm_preproc_write8.m

K   = size(res.mn,2);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3 = 1:d(3);

% prepare DCT parameters for bias correction
chan = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
d3      = [size(res.Tbias{1}) 1];
chan.B3 = spm_dctmtx(d(3),d3(3),x3);
chan.B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
chan.B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
chan.T  = res.Tbias{1};

% output image
Ysimu = zeros(d, 'single');

spm_progress_bar('init',length(x3),['Working on ' name],'Planes completed');
for z = 1:length(x3)

  % Bias corrected image
  f  = spm_sample_vol(res.image(1),x1,x2,o*x3(z),0);
  bf = exp(transf(chan.B1,chan.B2,chan.B3(z,:),chan.T));
  cr{1} = bf.*f;

  msk = any((f==0) | ~isfinite(f),3);

  % Parametric representation of intensity distributions
  q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
  q1  = reshape(q1,[d(1:2),numel(res.mg)]);
  
  % replace classes and means for GM/WM/CSF by external segmentations
  for k=1:3
    ind = find(res.lkp==k);
    for j=1:numel(ind)
      q1(:,:,ind(j)) = vol_seg(:,:,z,k);
      res.mn(1,ind(j)) = mn(k);
    end
  end

  % add WMHs as last class
  if ~isempty(WMH)
    q1(:,:,K) = 5*WMH(:,:,z);
  end
  
  s   = sum(q1,3);
  tmp = zeros(d(1:2));

  % sum up over first 3 classes and normalize to sum of 1
  for k=1:find(res.lkp<=3, 1, 'last')
    tmp = tmp + res.mg(k)*res.mn(1,k)*q1(:,:,k)./s;
  end

  % add remaining 3 BG classes from bias corrected image
  ind = tmp == 0;
  tmp(ind) = cr{1}(ind);
  tmp(msk) = 1e-3;

  Ysimu(:,:,z) = tmp;
  
  spm_progress_bar('set',z);
end
spm_progress_bar('clear');


%==========================================================================
% function [Ysimu, rf_field]  = add_bias_field(Ysimu, rf, idef_name, pth)
%
% Purpose
%   Apply a predefined RF bias field (A/B/C) from template space to native
%   space, scaled by rf.percent, optionally invert for negative percent.
%
% Inputs
%   Ysimu     - simulated image.
%   rf        - struct with .percent (signed), .type ('A'|'B'|'C').
%   idef_name - inverse deformation field for warping the RF template.
%   pth       - path to directory containing rf100_*.nii fields.
%
% Outputs
%   Ysimu     - modulated image.
%   rf_field  - applied RF field in native space.
%==========================================================================
function [Ysimu, rf_field] = add_bias_field(Ysimu, rf, idef_name, pth)

fprintf('Transform RF field to native space.\n');
% warp defined rf field to native space
rf_name = fullfile(pth,['rf100_' rf.type '.nii']);
rf_field = cat_vol_defs(struct('field1',{{idef_name}},'images',{{rf_name}},'interp',1,'modulate',0));
rf_field = single(rf_field{1}{1});

% apply defined percent and strength
rf_field = abs(rf.percent)/100 * (single(rf_field));

% invert field for neg. values
if rf.percent < 0
  rf_field = 1./rf_field;
end

ind = isfinite(rf_field);
rf_field = 1 + rf_field - mean(rf_field(ind));
Ysimu(ind) = rf_field(ind).*Ysimu(ind);


%==========================================================================
% function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf)
%
% Purpose
%   Create a smooth, random RF bias field via FFT-domain filtering with a
%   strength-dependent grid size, interpolate to image size, scale by rf.percent,
%   and apply to the simulated image.
%
% Inputs
%   Ysimu  - image to modulate.
%   rf     - struct: .type = [strength, rngSeed], .percent amplitude (signed).
%
% Outputs
%   Ysimu    - modulated image.
%   rf_field - generated RF field after smoothing and scaling.
%==========================================================================
function [Ysimu, rf_field] = add_simulated_bias_field(Ysimu, rf)

dim = size(Ysimu);

% set seeds to defined value
rng(rf.type(2),'twister')

N = 2^round(rf.type(1)); % Define the size of the 3D field w.r.t. defined strength
fwhm = 30; % Smoothing size
pad  = 3*fwhm; % pad border to prevent smoothing issues at image borders

% Generate a random 3D field
field = rand(N, N, N);

% Apply 3D FFT to the field
fieldFFT = fftn(field);

% Create a 3D frequency filter to manipulate smoothness
[x, y, z] = meshgrid(-N/2:N/2-1, -N/2:N/2-1, -N/2:N/2-1);
radius = sqrt(x.^2 + y.^2 + z.^2); % Distance from the center in 3D
smoothnessFilter = exp(-radius./(N/8)); % Gaussian filter for smoothness

% Shift the FFT for filtering
fieldFFTShifted = fftshift(fieldFFT);

% Apply the smoothness filter in the frequency domain
filteredFFT = fieldFFTShifted .* smoothnessFilter;

% Shift back and apply inverse 3D FFT
filteredField = ifftn(ifftshift(filteredFFT));

% Original grid
[x, y, z] = ndgrid(1:N, 1:N, 1:N);

% New grid dimensions
[xq, yq, zq] = ndgrid(linspace(1, N, dim(1)), ...
                      linspace(1, N, dim(2)), ...
                      linspace(1, N, dim(3)));

% Interpolate using interpn
filteredField = interpn(x, y, z, filteredField, xq, yq, zq, 'linear');

% resize field and scale it to 0..1
filteredField = filteredField - min(filteredField(:));
filteredField = filteredField/max(filteredField(:));

% pad border and create background 0f 0.5 to prevent smoothing issues at image border
vol = 0.5*ones(dim(1)+2*pad,dim(2)+2*pad,dim(3)+2*pad);
vol(pad+1:dim(1)+pad,pad+1:dim(2)+pad,pad+1:dim(3)+pad) = filteredField;

% smooth and resize back
spm_smooth(vol,vol,[fwhm fwhm fwhm]);
rf_field = vol(pad+1:dim(1)+pad,pad+1:dim(2)+pad,pad+1:dim(3)+pad);

% scale to an amplitude and mean of 1
rf_field = rf_field - min(rf_field(:));
rf_field = rf_field/max(rf_field(:));

% apply defined percent
rf_field = abs(rf.percent)/100 * rf_field;

% invert field for neg. values
if rf.percent < 0
  rf_field = 1./rf_field;
end
rf_field = 1 + rf_field - mean(rf_field(:));

% and finally apply bias field
Ysimu = rf_field.*Ysimu;


%==========================================================================
% function [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name)
%
% Purpose
%   Add simulated white matter hyperintensities (WMHs) consistent with a
%   prior WMH probability map and subject anatomy. The routine warps a WMH
%   atlas to native space, modulates it by a random field and a user-defined
%   strength, and injects the resulting WMH class into the PVE label image
%   and mixture model used for synthesis.
%
% Inputs
%   simu         - struct with field .WMH (double): strength parameter (>1
%                  increases WMH expression; values close to 1 keep it mild).
%   res          - SPM segmentation result struct (contains image header and
%                  GMM parameters mg/mn/vr/lkp; mn is updated for WMH class).
%   label_pve    - single(dims): current PVE-like label image in [1..3] (CSF=1,
%                  GM=2, WM=3). Will be extended to [1..4] where 4 encodes WMH.
%   template_dir - path containing the WMH prior map 'cat_wmh_miccai2017.nii'.
%   idef_name    - inverse deformation field to warp the WMH atlas to native
%                  space (continuous interpolation for intensities).
%
% Outputs
%   WMH          - single(dims): simulated WMH map in [0..1].
%   res          - struct: mixture model with an added WMH class (mn updated).
%   label_pve    - single(dims): updated label image where WMH contributes as
%                  an additional class (value 4 contribution), later used for
%                  synthesis.
%
% Algorithm
%   1) Warp the MICCAI2017 WMH prior map to native space and lightly smooth it.
%   2) Create a random 3D field, resample it to image dimensions, threshold to
%      form spatial support for patchy WMH distribution.
%   3) Erode WM to ensure spacing from GM and constrain WMHs to deep WM.
%   4) Combine WMH prior^(1/(strength-0.8)) with random field and WM mask,
%      smooth and normalize to [0,1].
%   5) Update label_pve by adding 2*WMH and clipping to max class label 4,
%      thereby introducing an additional WMH class contribution.
%   6) Extend the GMM by adding a WMH class intensity equal to the GM mean
%
% Notes
%   - Class encoding after this step: CSF=1, GM=2, WM=3, WMH=4.
%   - The strength parameter shapes the prior map nonlinearly; larger values
%     emphasize regions with higher WMH probability.
%   - The WM erosion step helps avoid spuriously labeling GM/CSF boundaries
%     as WMH.
%==========================================================================
function [WMH, res, label_pve] = simulate_WMHs(simu, res, label_pve, template_dir, idef_name)

Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
strength = simu.WMH;
vx = sqrt(sum(res.image(1).mat(1:3,1:3).^2));

% warp WMH map to native space
fprintf('Transform WMH map to native space.\n');
WMH_name = fullfile(template_dir,'cat_wmh_miccai2017.nii');
WMH = cat_vol_defs(struct('field1',{{idef_name}},'images',{{WMH_name}},'interp',1,'modulate',0));
WMH = single(WMH{1}{1});
dim = size(WMH);

% slightly smooth WMH atlas
spm_smooth(WMH, WMH, 2./vx);

% use filename for obtaining an integer that is used as RNG seed
[~, fname] = fileparts(res.image(1).fname);
seed = sum(double(fname));
rng(seed, 'twister')

N = 2^round(5); % Define the size of the initial 3D field

% Generate a random 3D field
fieldN = rand(N, N, N);

% Original grid
[x, y, z] = ndgrid(1:N, 1:N, 1:N);

% New grid dimensions
[xq, yq, zq] = ndgrid(linspace(1, N, dim(1)), ...
                      linspace(1, N, dim(2)), ...
                      linspace(1, N, dim(3)));

% Interpolate using interpn
field = interpn(x, y, z, fieldN, xq, yq, zq, 'linear');
field = single(field>0.7);

% find WM in original image and erode it to ensure some space to GM
WM = cat_vol_morph(Yp0toC(label_pve, 3) > 0.5,'de',3,vx);

% apply strength parameter to WM atlas and combine it with field and WM
WMH = WMH.^(1/(strength-0.8)).*field.^2.*WM;
spm_smooth(WMH, WMH, 2./vx);
WMH = WMH/max(WMH(:));

% correct PVE label and add additional WMH class (which should be a
% increased in intensities)
label_pve = label_pve + 5*WMH/strength.^0.75;
label_pve(label_pve > 4) = 4;

% use mean of GM for additional class
K   = numel(res.mg);
intensity_WMH = mean(res.mn(1,res.lkp==1));
K = K + 1;
res.mn(1,K) = intensity_WMH;


% function Yb = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th)
% SKULL_STRIP_APRG - Brain extraction via CAT12 APRG
%
% Purpose
%   Compute a robust brain mask using CAT12's Adaptive Probability Region-Growing
%   (APRG) method, leveraging SPM/CAT tissue posteriors and data-driven
%   intensity thresholds from the current image. This improves downstream PVE
%   label generation by restricting to brain voxels only.
%
% Inputs
%   Ysrc - single(dims): Source image used for threshold estimation (same space
%          as SPM/CAT outputs from preproc). Typically the first channel.
%   Ycls - 1x6 cell of uint8(dims): Tissue class posteriors (0..255) from
%          cat_spm_preproc_write8, in SPM/CAT convention (c1..c6).
%   res  - struct: SPM/CAT segmentation result with fields used here:
%              .mn  (means per class/component)
%              .mg  (mixture weights)
%              .lkp (lookup of class index per component, 1..6)
%   dim  - [nx ny nz]: Volume dimensions of the images.
%   T3th - 1x3 double: Intensity thresholds [CMth, CSFth, WMth] used by APRG
%          (derived from get_tissue_thresholds). CMth is a central-matter
%          threshold between GM and WM, CSFth approximates CSF mean, WMth is a
%          robust WM anchor.
%
% Output
%   Yb   - logical/single(dims): Brain mask estimated by APRG (1=brain, 0=non-brain).
%
% Algorithm (summary)
%   1) Estimate class-specific intensity anchors via weighted means of the GMM
%      parameters (clsint for CSF/GM/WM).
%   2) Build a 4D stack P(:,:,:,1..6) from Ycls (uint8 posteriors), as required
%      by CAT12’s APRG routine.
%   3) Derive robust thresholds for WM and a central matter (CM) value using the
%      WM posterior median within Ysrc, which mitigates issues in presence of WMHs.
%   4) Call cat_main_APRG(Ysrc, P, res, T3th) to obtain the brain mask.
%
% Notes
%   - Assumes T1-like contrast (WM > GM > CSF) for threshold reasoning; a median
%     WM intensity is used to be robust to WMHs.
%   - Requires CAT12 on the MATLAB path; uses cat_main_APRG.
%   - The output mask matches the input volume dimensions and orientation.
%==========================================================================
function Yb = skull_strip_APRG(Ysrc, Ycls, res, dim, T3th)

% we need a 4D array for cat_main_APRG
P = zeros([dim 6],'uint8');
for i=1:6
  P(:,:,:,i) = Ycls{i};
end

res.isMP2RAGE = 0;
Yb = cat_main_APRG(Ysrc, P, res, T3th);


function T3th = get_tissue_thresholds(Ysrc, Ycls, res)
% GET_TISSUE_THRESHOLDS - Estimate robust CSF/GM/WM thresholds
%
% Purpose
%   Compute three intensity anchors used by LAS/APRG steps:
%     - WMth: a robust white-matter threshold derived from the weighted WM
%       mean (from the GMM) and the median intensity within high WM posterior
%       voxels to mitigate WMH effects.
%     - CSFth: an approximation of the CSF intensity anchor from the GMM.
%     - CMth: a central-matter threshold between GM and WM that adapts to the
%       current image contrast.
%
% Inputs
%   Ysrc - single(dims): Source image used to compute medians within WM.
%   Ycls - 1x6 cell of uint8(dims): SPM/CAT posteriors (0..255), where
%          Ycls{2} corresponds to WM.
%   res  - struct: SPM/CAT segmentation result providing Gaussian mixture
%          parameters mn, mg, and lookup lkp.
%
% Output
%   T3th - 1x3 double: [CMth, CSFth, WMth] thresholds used by skull stripping
%          and LAS normalization routines.
%
% Algorithm
%   1) clsint(k): weighted class mean from the GMM for class k.
%   2) WMth: max(clsint(WM), median(Ysrc within high WM posterior)).
%   3) CMth: if intensities are inverted (CSF > WM), use CSF anchor; otherwise
%      use 2*GMmean - WMth, clipped not to exceed CSF anchor.
%   4) Return [CMth, CSFmean, WMth].
%
% Notes
%   - Class indices follow SPM/CAT convention (typically: 1=GM, 2=WM, 3=CSF).
%   - High WM posterior is defined with a conservative threshold (Ycls{2}>192
%     on 0..255 scale) to obtain a stable median.
%   - These anchors are designed for T1-like contrast but include a simple
%     inversion check (CSF > WM) for robustness.

clsint = @(x) round( sum(res.mn(res.lkp==x) .* res.mg(res.lkp==x)') * 10^5)/10^5;

% Use median for WM threshold estimation to avoid problems in case of WMHs
WMth = double(max(clsint(2), cat_stat_nanmedian(Ysrc(Ycls{2}>192)))); 
if clsint(3)>clsint(2) % invers
  CMth = clsint(3); 
else
  CMth = min([clsint(1) - diff([clsint(1),WMth]), clsint(3)]);
end

T3th = double([CMth, clsint(1), WMth]);