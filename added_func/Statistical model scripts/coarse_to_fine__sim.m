%% what happens for fixed  Nscales Ntheta for different scales shifts
%% test model script
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd (myPath);
cd ('../../');
initpath;
cd(myPath);
%% create data:
load volume.mat
n_theta_vec=[360 180 120 90];
max_scales = linspace(1.12,2,6) ;
max_shifts = linspace(5,15,6);

[Nx, Ny, Nz] = size(vol);
n=5;
inv_rot_matrices = zeros(3, 3, n);
ref_q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
for k = 1:n
    quat = ref_q(:, k);
    ref_q(:, k) = quat/norm(quat);
    rot = q_to_rot(ref_q(:, k));
    inv_rot_matrices(:, :, k) = rot';
end
projections=cryo_project(vol,q_to_rot(ref_q),size(vol,1),'single'); % generate phantom projecitons

% The following step is ESSENTIAL before calling FIRM reconstruction on the
% projections.
projections=permute(projections,[2 1 3]);   % transpose each image

% noise_vec=[1 100];
clstack=zeros(n,n,length(max_shifts),length(max_scales),length(n_theta_vec));

ref_clmatrix=clmatrix_cheat(q_to_rot(ref_q),360);
for ii = 1:length(max_shifts)
    for j = 1:length(max_scales)
        for kk= 1:length(n_theta_vec)
            n_theta = n_theta_vec(kk);
            max_scale = max_scales(j);
            max_shift = max_shifts(ii);
            
            % Pad projections
            pad_factor = 1.15*max_scale;
            
            imsize = size(projections,1);
            padded_size = floor(pad_factor*[imsize, imsize]);
            
            padded_projections = padarray(projections, ceil((padded_size-imsize)/2), 'pre');
            padded_projections = padarray(padded_projections, ceil((padded_size-imsize)/2), 'post');
            
            max_log_scale = log(max_scale);
            
            log_scales = 2*max_log_scale*rand(n,1) - max_log_scale;
            log_scales = log_scales - mean(log_scales);
            scales = exp(log_scales);
            
            scaled_projections = scale_projections(padded_projections, scales);
            
            % Add shifts to projections
            [scaled_projections, ref_shifts] = cryo_addshifts(scaled_projections,[],max_shift,1);
            %projections=scaled_projections; ref_shifts=0;
            scaled_projections(scaled_projections<0)=0;
            max_proj = max(scaled_projections(:));
            max_well = 20000;
            dark_current = 10 * 1e-6 * ones(size(scaled_projections(:,:,1)));
            % Add noise
            for i = 1:n
                scaled_projections(:,:,i) = max_proj/(max_well*0.5*1e-6) * ...
                    (imnoise(single(dark_current),'poisson') + imnoise(imnoise(0.5*max_well*1e-6*single(scaled_projections(:,:,i)/max_proj), 'poisson'),'speckle', 1e-6));
            end
            
            fprintf('Commonline detection for noisy shifted and scaled images using Gaussian filter.\n');
            fprintf('=====================================================================\n');
            mask_radius = round(0.75*size(scaled_projections,1));
            n_r = round(0.75*size(scaled_projections,1));
            Nscales=10;
            Nshifts=10;
            [masked_projs,~]=mask_fuzzy(scaled_projections,mask_radius);        %% Mask projections
            log_message('Computeing polar Fourier transform of projections. n_theta=%d, n_r=%d',n_theta,n_r);
            [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
            [clstack(:,:,ii,j,kk)] = only_commonlines_gaussian_shift_scale(pf,max_shift,Nshifts, Nscales, max_scale);%% Find only common lines from projections
            save('coarse_to_fine__sim');
        end
    end
end

