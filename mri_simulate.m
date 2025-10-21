function mri_simulate(simu, rf)
% MRI_SIMULATE - Simulates MR images with various optional features
%
% Overview:
%   `mri_simulate` generates simulated MRI images, allowing for T1-weighted
%   (T1w) imaging simulations. It employs a high-resolution
%   (e.g., 0.5mm) T1w image, typically Colin27, as a base. Users can introduce
%   various artifacts and features like Gaussian noise, and RF B1 inhomogeneities. 
%   It supports simulations of atrophy or cortical thickness modifications. 
%   Preprocessing with SPM12 segmentation is required for custom images.
%
% Syntax:
%   mri_simulate(simu, rf)
%
% Parameters:
%   simu (struct): Simulation parameters.
%       - 'name' (string): Filename of the T1w input image.
%       - 'pn' (double): Percentage noise level to introduce Gaussian noise.
%       - 'rng' (double or []): Seed for the random number generator; use []
%         for MATLAB's default behavior.
%       - 'resolution' (double or [x, y, z]): Spatial resolution of the
%         simulated image. Use 'NaN' for keeping original image resolution
%       - 'atrophy' (cell): Specifies regions of interest (ROIs) for simulating
%         atrophy, including atlas name, ROI IDs, and atrophy values. Multiple ROIs 
%         and the respective atrophy values can be defined. An atrophy value of 
%         1.5 leads to a GM reduction of about 5%, while a value of 2 corresponds 
%         to around 10% and 3 to around 15%. Please note, that this option takes 
%         a lot of time because the atlas labels have to be interpolated using 
%         categorical interpolation (i.e. each label separetely). Either thickness 
%         or atrophy can be simulated.
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
%       - 'save' (logical): Flag to save the ground truth label if set to 1.
%   rf (struct): RF bias field parameters.
%       - 'percent' (double): Amplitude of the bias field in percentage.
%         Negative values invert the field.
%       - 'type' (char or [int, int]): Specifies the bias field type, options are 
%         'A', 'B', or 'C' for predefined fields from MNI or an integer array with 
%         two numbers. The first integer value adjusts the local strength of the 
%         field by varying maximum frequencies of the FFT. Meaningful values are 
%         1..4, while a value of 3 or 4 corresponds to a bias field of a 7T scanner 
%         (without further correction such as in mp2rage). The second integer sets 
%         the random generator to that seed value, which allows simulating different 
%         bias fields.
%       - 'save' (logical): Option to save the simulated bias field if 'type' is
%         numeric and 'save' is set to 1.       
%
% Optional Inputs:
%   Default values are used if 'simu' or 'rf' parameters are not provided. See
%   examples for default structures.
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
%   Example 1 - Basic simulation with specific noise:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3,...
%                     'resolution', 0.5,...
%                     'atrophy', [], 'rng', 42);
%       rf = struct('percent', 20, 'type', 'A','save',0);
%       mri_simulate(simu, rf);
%
%   Example 2 - Advanced simulation with atrophy and custom RF field and 
%               large slice thickness:
%       simu = struct('name', 'custom_t1.nii', 'pn', 2,...
%                     'resolution', [0.5, 0.5, 1.5],...
%                     'rng', []);
%       simu.atrophy = {'hammers',[28, 29], [2, 3]};
%       rf = struct('percent', 15, 'type', [3, 42]);
%       mri_simulate(simu, rf);
%
%   Example 3 - Thickness simulation:
%       simu = struct('name', 'colin27_t1_tal_hires.nii', 'pn', 3,...
%                     'resolution', NaN,...
%                     'atrophy', [], 'rng', [],...
%                     'thickness', [1.5 2.0 2.5]);
%       rf = struct('percent', 20, 'type', 'A');
%       mri_simulate(simu, rf);
%
% TODO: simulation of motion artefacts using FFT and shift of phase information

% Default simulation parameters
def.name       = 'colin27_t1_tal_hires.nii';
def.pn         = 3;
def.resolution = NaN;
def.atrophy    = [];
def.thickness  = 0;
def.rng        = 0;
def.save       = 1;

if nargin < 1, simu = def;
else, simu = cat_io_checkinopt(simu, def); end

% Default bias field parameters
def.percent    = 20;
def.type       = [2 0];
def.save       = 0;

if nargin < 2, rf = def;
else, rf = cat_io_checkinopt(rf, def); end

[pth, name, ~] = spm_fileparts(simu.name);
pth_root = fileparts(which(mfilename));

% name of seg8.mat file that contains SPM12 segmentation parameters
mat_name = fullfile(pth, [name '_seg8.mat']);

% CAT12 template dir is later used for defining atrophy atlas
template_dir = fullfile(spm('dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym');

% call SPM segmentation if necessary
if ~exist(mat_name,'file')
  fprintf('We have to run SPM segmentation first.\n')

  matlabbatch{1}.spm.spatial.preproc.channel.vols = {simu.name};
  spm_jobman('run',matlabbatch);
  clear matlabbatch
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

% load seg8.mat file and define some parameters
res = load(mat_name);
[~, bname, ext] = spm_fileparts(res.image(1).fname);
res.image(1) = spm_vol(fullfile(pth,[bname ext]));
V = res.image(1);
d   = V.dim(1:3);
vx = sqrt(sum(V.mat(1:3,1:3).^2));

% obtain SPM segmentations
[Ysrc, Ycls, Ydef] = cat_spm_preproc_write8(res,zeros(max(res.lkp),4),zeros(1,2),[0 0],0,2);

% write forward deformation field
idef_name = fullfile(pth, ['y_inv_' name ext]);

if ~exist(idef_name,'file')
  def_name = fullfile(pth, ['y_' name ext]);
  Ndef = nifti;
  Ndef.dat  = file_array(def_name,[d,1,3],...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
  Ndef.mat = V.mat;
  Ndef.mat0 = V.mat;
  create(Ndef);
  Ndef.dat(:,:,:,1,1) = Ydef(:,:,:,1);
  Ndef.dat(:,:,:,1,2) = Ydef(:,:,:,2);
  Ndef.dat(:,:,:,1,3) = Ydef(:,:,:,3);
  
  % write inverse deformation field
  idef_name = invert_deformations(simu);
  spm_unlink(def_name);
  clear Ydef
end

% LAS correction and SANLM denoising
Ycorr = cat_main_LASsimple(Ysrc,Ycls);
cat_sanlm(Ycorr,3,1);

% Replace GM/WM/CSF segmentation by labels using LAS corrected image
Yp0toC = @(Yp0,c) 1-min(1,abs(Yp0-c));
Yseg = zeros([d, 3]);
brainmask = zeros(d);
seg_order = [2 3 1];
for i = 1:3
  brainmask = brainmask + double(Ycls{i})/255;
end
for i = 1:3
    Yseg(:,:,:,i) = Yp0toC(3*Ycorr, seg_order(i)).*brainmask; 
end
clear Ycls

% add atrophy to GM by decreasing GM value in ROI and increasing CSF value
if simu_atrophy
  Yseg = simulate_atrophy(simu, Yseg, d, template_dir, idef_name);
end

% get ground truth label using GM/WM/CSF
label_pve = zeros(d, 'single');
order = [3 1 2];
for k = 1:3
    label_pve = label_pve + k*(Yseg(:,:,:,order(k)));
end

% get means for GM/WM/CSF
mn = zeros(3,1);
for k=1:3
  ind = find(res.lkp==k);
  for j=1:numel(ind)
    mn(k) = mn(k) + res.mg(ind(j))*res.mn(1,ind(j));
  end
end

% extend target voxel size if needed
if isscalar(simu.resolution)
  if isnan(simu.resolution)
    simu.resolution = vx;
  else
    simu.resolution = simu.resolution * ones(1,3);
  end
else
  if isnan(simu.resolution(1))
    simu.resolution = vx;
  end
end

if any(simu.thickness)
  [label_pve, Yseg] = simulate_thickness(label_pve, simu, Yseg, d, template_dir, idef_name, vx, order);
end

Ysimu = synthesize_from_segmentation(Yseg, name, res, mn, d);

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
relWM_peak = mn(2)/mx_vol;

% output matrix
Vres.dim = round(V.dim.*vx./simu.resolution);
P = spm_imatrix(V.mat);
P(7:9) = P(7:9)./vx.*simu.resolution;
P(1:3) = P(1:3) + vx - simu.resolution;
Vres.mat = spm_matrix(P);

% output in defined resolution
volres   = zeros(Vres.dim);
if simu.save
  labelres_pve = zeros(Vres.dim);
end

for sl = 1:Vres.dim(3)
  M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
  M1 = Vres.mat\V.mat\M;
  
  % use sinc interpolation for simulated image
  volres(:,:,sl) = spm_slice_vol(Ysimu,M1,Vres.dim(1:2),-5);
  if simu.save
    % and linear interpolation for label image
    labelres_pve(:,:,sl) = spm_slice_vol(label_pve,M1,Vres.dim(1:2),1);
  end
end

volres = volres / mx_vol;
% get higher peak of WM for estimating noise level
noise_std_dev = relWM_peak * simu.pn / 100;

% optionally ensure same noise for every trial
if ~isempty(simu.rng)
  rng(simu.rng,'twister');
end

% generate and add Gaussian noise and correct potential neg. values
noise = (noise_std_dev) * randn(size(volres));
volres = volres + noise;
volres = mx_vol * max(min(volres, 1), 0);

str = '';
if (exist('rf','var') && rf.percent ~= 0)
  if isnumeric(rf.type)
    str = sprintf('_rf%g_%d_%d', rf.percent, rf.type(1), rf.type(2));
  else
    str = sprintf('_rf%g_%s', rf.percent, rf.type);
  end
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
mind = sum(Yseg(:,:,:,1:3),4) > 0.15;
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
spm_write_vol(Vres, volres.*mind);

% write ground truth label
if simu.save
  if any(simu.thickness)
    label_pve_name = fullfile(pth, sprintf('label_pve_%gmm_%s%s.nii', mean_resolution, name, str));
  else
    label_pve_name = fullfile(pth, sprintf('label_pve_%gmm_%s.nii', mean_resolution, name));
  end
  
  fprintf('Save %s\n', label_pve_name);
  Vres.fname = label_pve_name;
  Vres.pinfo = [1/255/3 0 352]';
  Vres.dt    = [4 0];
  spm_write_vol(Vres, labelres_pve);

end

% save simulated bias field if defined
if rf.save
  rf_name = fullfile(pth, sprintf('%s_%gmm_%s.nii', str(2:end), mean_resolution, name));
  fprintf('Save %s\n', rf_name);
  Vres.fname = rf_name;
  Vres.pinfo = [1/max(rf_field(:)) 0 352]';
  Vres.dt    = [2 0];
  spm_write_vol(Vres, rf_field);
end


%==========================================================================
% function t = transf(B1,B2,B3,T)
% from spm_preproc_write8.m
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
% function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, order)
%==========================================================================
function [label, Yseg] = simulate_thickness(label, simu, Yseg, d, template_dir, idef_name, vx, order)

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
label0 = label;

label1 = cell(numel(simu.thickness),1);

Yseg(:,:,:,1:3) = 0;

% vary range of PVE from -0.25..0.25 in 20 steps
pve_range = linspace(-0.25,0.25,20);
for pve_step = 1:numel(pve_range)

  % euclidean distance to wm
  wm  = round(label+pve_range(pve_step)) > gwm_val;
  wm = cat_vol_morph(wm,'close',1);
  
  % euclidean distance to wm
  Dwm = (bwdist(wm) - 0.5) * mean(vx); % also consider voxelsize/2 correction of distance

  for k=1:numel(simu.thickness)
  
    label1{k} = round(label+pve_range(pve_step));
  
    % set gm and csf to csf
    label1{k}(label1{k} < wm_val & label1{k} > csf_val) = csf_val;
  
    % set dilated gm to gm
    label1{k}(label1{k} < wm_val & Dwm <= simu.thickness(k)) = gm_val;
      
    % check distance inside gm
    D = (bwdist(label1{k}==csf_val) + bwdist(label1{k}==wm_val) - 1)*mean(vx);
    D(label1{k}~=gm_val) = 0;
        
    % dilate gm where distance in gm is too large to make more space
    gm1 = cat_vol_morph(D > 0.97*simu.thickness(k),'o',1);
    gm1 = cat_vol_morph(gm1,'dd',2);
    gm1(~mask) = 0;
    label1{k}(gm1 > 0 | label1{k} == gm_val) = csf_val;
    
    % euclidean distance to wm
    wm  = label1{k}==wm_val;
    Dwm = (bwdist(wm) - 0.5) * mean(vx);
    
    label1{k}(~wm) = csf_val;
    label1{k}(~mask) = 0;
    
    % set dilated gm to gm
    label1{k}(label1{k} == csf_val & Dwm <= simu.thickness(k)) = gm_val;
  
    % get original label for basal ganglia and cerebellum
    label1{k}(mask_orig) = label0(mask_orig);
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
% function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d)
%==========================================================================
function Ysimu = synthesize_from_segmentation(vol_seg, name, res, mn, d)
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
  
  s   = sum(q1,3);
  tmp = zeros(d(1:2));

  % sum up over all classes and normalize to sum of 1
  for k=1:K
    tmp = tmp + res.mn(1,k)*q1(:,:,k)./s;
  end
  tmp(msk) = 1e-3;
  Ysimu(:,:,z) = tmp;
  
  spm_progress_bar('set',z);
end
spm_progress_bar('clear');


%==========================================================================
% function [Ysimu, rf_field]  = add_bias_field(Ysimu, rf, idef_name, pth)
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
% function idef_name = invert_deformations(simu)
%==========================================================================
function idef_name = invert_deformations(simu)

% run deformation utility if inverse deformation field does not exist
[pth, name, ext] = spm_fileparts(simu.name);
if ~exist(fullfile(pth,['y_inv_' name ext]),'file')
  matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {fullfile(pth,['y_' name ext])};
  matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {simu.name};
  matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = ['inv_' name];
  matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {pwd};
  spm_jobman('run',matlabbatch);
  clear matlabbatch
end
idef_name = fullfile(pth,['y_inv_' name ext]);
