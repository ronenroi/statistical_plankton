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
cd('..\..\');
initpath;
cd(myPath);

%% load volume
load Bgal_3d.mat
vol = Bgal_3d;
%% Stretch projections by a random factor 
% max_scales = [1.0, 1.15, 1.25, 1.5, 1.65, 1.75, 2, 2.25, 2.5, 2.75, 3];
max_scales = 2;
[Nx, Ny, Nz] = size(vol);

for j = 1:length(max_scales)
    max_scale = max_scales(j);
    % Generate random rotations.
    % Create 100 random uniformly sampled rotations (from the uniform
    % distribution on SO(3)).

    n=100;                              % use 500 images
    inv_rot_matrices = zeros(3, 3, n);
    q = qrand(n);                       % generate rotations as quaternions
    % find inverse rotation matrices
    for k = 1:n
        quat = q(:, k);
        q(:, k) = quat/norm(quat);
        rot = q_to_rot(q(:, k));
        inv_rot_matrices(:, :, k) = rot';
    end

    %% Generate 2D projections.
    % Caclulate the 2D proejctions of the 3D volume in the directions
    % determined by the randomly generated rotations matrices.

    projections=cryo_project(vol,q_to_rot(q),size(vol,1),'single'); % generate phantom projecitons

    % The following step is ESSENTIAL before calling FIRM reconstruction on the
    % projections.
    projections=permute(projections,[2 1 3]);   % transpose each image

    % Pad projections 
    pad_factor = 1.15*max_scale;
    
    imsize = size(projections,1);
    padded_size = floor(pad_factor*[imsize, imsize]);

    padded_projections = padarray(projections, ceil((padded_size-imsize)/2), 'pre');
    padded_projections = padarray(padded_projections, ceil((padded_size-imsize)/2), 'post');

    scales = (max_scale - 1/max_scale)*rand(n,1) + 1/max_scale;
    mu_scales = mean(scales);
    [X,Y,Z] = meshgrid(linspace(1,Nx,mu_scales*Nx), ...
                       linspace(1,Ny,mu_scales*Ny), ...
                       linspace(1,Nz,mu_scales*Nz));
                   
    vol_upscaled = interp3(vol,X,Y,Z);
    scaled_projections = scale_projections(padded_projections, scales);

    % Add shifts to projections
    max_shift = 10;
    [projections, ref_shifts] = cryo_addshifts(scaled_projections,[],max_shift,1);
    %projections=scaled_projections; ref_shifts=0;
    projections(projections<0)=0;
    max_proj = max(projections(:));
    max_well = 20000;
    dark_current = 10 * 1e-6 * ones(size(projections(:,:,1)));
    % Add noise 
    for i = 1:n
        projections(:,:,i) = max_proj/(max_well*0.5*1e-6) * ...
            (imnoise(single(dark_current),'poisson') + imnoise(imnoise(0.5*max_well*1e-6*single(projections(:,:,i)/max_proj), 'poisson'),'speckle', 1e-6));
    end


    % Save all necessary matrices
    ref_q = q;
    ref_inv_rot_matrices = inv_rot_matrices;
    ref_scales = scales;
    ref_vol = vol_upscaled;
    fprintf(['Saving matrices ' num2str(j) '/' num2str(length(max_scales)) '\n']);
    save(['simulated_data\scale' num2str(max_scale) '_Bgal_3d.mat'], 'ref_q', 'ref_inv_rot_matrices',...
        'projections', 'ref_scales', 'ref_vol', 'ref_shifts', 'max_scale', 'max_shift') 
end