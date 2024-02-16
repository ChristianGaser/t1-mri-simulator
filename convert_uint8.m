function convert_uint8

P = spm_select(Inf,'image','select images');
V = spm_vol(P);
flags = struct('dmtx', 0, 'mask', 0, 'interp', 1, 'dtype', 2);

for i = 1:numel(V)
  spm_imcalc(V(i),'output.nii','i1',flags);
  system(['mv output.nii ' V(i).fname]);
end