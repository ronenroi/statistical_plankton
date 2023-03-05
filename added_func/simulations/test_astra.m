%% Using simulated projections
%
% Usage example of functions for generating simulated projections.
%
% The script generates a volume sampled from a 3D ribosome model,
% computes 2D projections of the sampled volume with
% differente translation, scale, rotation and noise
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd (myPath);
cd('..\..\');
initpath;
cd('..\');
astradir = cd;
addpath(genpath(fullfile([astradir,'\astra-1.9.0.dev11'])));
cd(myPath);

%% load volume
load volume.mat
vol(vol<0)=0;
%% Stretch projections by a random factor
% max_scales = [1.0, 1.15, 1.25, 1.5, 1.65, 1.75, 2, 2.25, 2.5, 2.75, 3];
[Nx, Ny, Nz] = size(vol);

% Generate random rotations.
% Create 100 random uniformly sampled rotations (from the uniform
% distribution on SO(3)).

n=100;                              % use 100 images
rot_matrices = zeros(3, 3, n);
inv_rot_matrices = zeros(3, 3, n);

q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
for k = 1:n
    quat = q(:, k);
    q(:, k) = quat/norm(quat);
    rot = q_to_rot(q(:, k));
    rot_matrices(:, :, k) = rot;
    inv_rot_matrices(:, :, k) = rot';
end


%% Generate 2D projections.
% Caclulate the 2D proejctions of the 3D volume in the directions
% determined by the randomly generated rotations matrices.

projections=cryo_project(vol,rot_matrices,size(vol,1),'single'); % generate phantom projecitons
projections=permute(projections,[1 3 2]);   % col x n_proj x row

%%
centered_projections = cryo_addshifts(projections,-gt_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./gt_scales);
projsCn=permute(unscaled_centered_projections,[2 3 1]);   % col x n_proj x row

vectors = rot2astraVec(inv_rotations);
vol_geom = astra_create_vol_geom(size(permute(gt_vol,[2,3,1])));
proj_geom = astra_create_proj_geom('parallel3d_vec', size(projsCn,3), size(projsCn,1), vectors);
%[proj_id, proj_data] = astra_create_sino3d_cuda(gt_vol, proj_geom, vol_geom);

% figure;viewstack(permute(projsCn,[1,3,2]),5,5)
% figure;viewstack(permute(proj_data,[1,3,2]),5,5)

[vol_id, recon_vol] = astra_create_backprojection3d_cuda(projsCn, proj_geom, vol_geom);
astra_mex_data3d('delete', vol_id);

%%
vectors = rot2astraVec(inv_rot_matrices);
vol_geom = astra_create_vol_geom(size(permute(vol,[2,3,1])));
proj_geom = astra_create_proj_geom('parallel3d_vec', size(projections,3), size(projections,1), vectors);

% Create projection data from this
[proj_id, proj_data] = astra_create_sino3d_cuda(double(vol), proj_geom, vol_geom);
figure;viewstack(permute(projections,[1,3,2]),5,5)
figure;viewstack(permute(proj_data,[1,3,2]),5,5)
rec_id = astra_mex_data3d('create', '-vol', vol_geom);
cfg = astra_struct('SIRT3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id, 150);
% Get the result
recon_vol = astra_mex_data3d('get', rec_id);

% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);
% %%
% [bproj_id, recon_vol] = astra_create_backprojection3d_cuda(projections, proj_geom, vol_geom);
% 
% astra_mex_data3d('delete', bproj_id);
%recon_vol = recon3D_astra(projections,inv_rot_matrices, size(projections,1));

%% 3D plots

vol_vec{1} = vol;
vol_vec{1}(vol<0)=0;
vol_vec{2} = recon_vol;
vol_vec{2}(recon_vol<0) = 0;


titles_vec{1} = 'Ground truth';
titles_vec{2} = 'Reconstruction with gt scale-shift';


figure;
for i=1:length(vol_vec)
    curr_vol = vol_vec{i};
    thr = max(vol_vec{i}(:)) * 10 / 100;
    err = norm(curr_vol(:)-vol_vec{1}(:))/norm(vol_vec{1}(:));
    curr_vol(curr_vol < thr) = 0;
    subplot(2,length(vol_vec)/2,i);
    vol3d('cdata',curr_vol);
    title({titles_vec{i},['relative error = ' num2str(err)]});
    colormap(gray(256));
    caxis([0 max(vol_vec{i}(:))])
    % alphamap(alpha);
    axis equal off
    set(gcf, 'color', 'w');
    view(3);
end

%% projections comparison
n=5;                              % projection number
q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
lim_max = -inf;
lim_min = inf;

for k = 1:n
    quat = q(:, k);
    q(:, k) = quat / norm(quat);
end
for i=1:length(vol_vec)
    curr_vol = vol_vec{i};
    proj = cryo_project(curr_vol,q_to_rot(q),size(curr_vol,1),'single'); % generate projecitons
    test_projections{i} = permute(proj,[2 1 3]);   % transpose each image
    lim_max = max(max(test_projections{i}(:)),lim_max);
    lim_min = min(min(test_projections{i}(:)),lim_min);
    
end
%% show projections
figure
colormap(gray)
for i=1:length(test_projections)
    curr_proj = test_projections{i};
    for j=1:size(curr_proj,3)
        subplot(length(test_projections),size(curr_proj,3),(i-1)*size(curr_proj,3)+j)
        curr_im = curr_proj(:,:,j);
        imagesc(curr_im);
        axis off
        axis image
        title({titles_vec{i},['proj. #' int2str(j)]});
    end
    caxis([lim_min lim_max])
    
end





