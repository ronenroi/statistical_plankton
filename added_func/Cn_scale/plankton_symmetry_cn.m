clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd (myPath);
cd ('../../');
initpath;
cd(myPath);
%%
if 1
recon_mat_fname = 'results/sym_vol_est';
    n_Images = 50;
    max_scale = 1.2;
if 1
[gt_vol,projections,inv_rotations,gt_q,gt_shifts,gt_scales,n_symm] = generate_symetric_projections(n_Images,max_scale);
save(recon_mat_fname,'gt_vol','projections','inv_rotations','gt_q','gt_shifts','gt_scales','n_symm');
else
    load(recon_mat_fname,'gt_vol','projections','inv_rotations','gt_q','gt_shifts','gt_scales','n_symm');
end
else
    recon_mat_fname = 'results/sym_plankton_vol_est_06_11_real';
    dataFolder = '../../../data sets/jules images/star/processed/usable'; projectionType = 'DarkField';  
 projections = func_loadPlanktonData(dataFolder, projectionType);
n_Images = size(projections,3);
n_symm = 8;
max_scale = 1.2;
save(recon_mat_fname);
end
n_theta = 360;
n_r = round(0.75*size(projections,1));
max_shift = 15;
shift_step = 0.5;
Nscales = 15;
mask_radius = round(0.75*size(projections,1));
inplane_rot_res = 0.5;
n_Points_sphere = 400;

cn_symmetry_vol_est_iter(recon_mat_fname,projections,n_symm,...
    n_theta,n_r,max_shift,shift_step,max_scale,Nscales,mask_radius,inplane_rot_res,n_Points_sphere);%,gt_q,gt_vol


