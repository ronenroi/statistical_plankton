%% Basic abinitio reconstruction - example 1.
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
register=1;
%% Create volume
load volume
%%
x0 = round(size(vol,1)/2);
if 0
sym_vol = zeros(size(vol));
sym_vol(x0:end,x0:end,x0:end) = vol(x0:end,x0:end,x0:end) ;
sym_vol = sym_vol + fliplr(sym_vol) ;
sym_vol = sym_vol + flipud(sym_vol) ;
sym_vol = sym_vol + flip(sym_vol,3) ;
figure;vol3d('Cdata',sym_vol);view(3);
end
flat_vol = vol(:,:,x0-2:x0+2);
%% create projections
vol = flat_vol;
max_scale = 2;
n=100;                              % use 100 images
    inv_rot_matrices = zeros(3, 3, n);
    q = qrand(n);                       % generate rotations as quaternions
    % find inverse rotation matrices
    for k = 1:n
        quat = q(:, k);
        q(:, k) = quat/norm(quat);
        rot = q_to_rot(q(:, k));
        inv_rot_matrices(:, :, k) = rot';
    end
        projections=cryo_project(vol,q_to_rot (q),size(vol,1),'single'); % generate phantom projecitons

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

viewstack(projections,5,5);   % Display the proejctions.

%% Compute common lines from projections

n_theta=180; % Angular resolution - number of sinograms computed for each 
            % projection. This corresponds to a resolution of 2 degrees.

% Find common lines between projections. Note that the max_shift used is
% much larger than above. This is because it is a 1D shift used to seearch
% for common lines. If d is the maximal expected 2D shift in a projection
% (in any direction), then the 1D shift used to find common lines should
% equal ceil(2*sqrt(2)*d).  
max_shift = ceil(2*sqrt(2)*max_shift) ;
shift_step = 1;
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));

[np,~]=mask_fuzzy(projections, mask_radius);
[npf,~]=cryo_pft(np,n_r,n_theta);
% [ clstack_SPR,~, shift_equations,~] = ...
%         commonlines_gaussian(npf,max_shift,shift_step ); %w.o scales and shifts estimation

% estimate shifts and scales with iterative refinement
clusteringThresh = 0.025;


[unscaled_centered_projections, iterative_LS_est_scales, ...
 ~, clstack_GSPR, ~, iterative_LS_est_shifts, ref_q_GSPR,valid_projections,valid_scales] = estimate_scales_shifts(projections, max_shift,...
                                                      max_scale, mask_radius, n_r, n_theta, clusteringThresh, q,scales);
%[est_inv_rots_LS_GSPR] = est_orientations_LS(clstack_GSPR, n_theta);
% [est_inv_rots_LUD_SPR] = est_orientations_LUD(clstack_SPR, n_theta); % Least Unsquared Deviations
[est_inv_rots_LUD_GSPR] = est_orientations_LUD(clstack_GSPR, n_theta); % Least Unsquared Deviations

[MSE_LUD_GSPR, ~, ~, aligned_rots_GSPR_LUD] = check_MSE(est_inv_rots_LUD_GSPR, ref_q_GSPR);
disp(['MSE LUD GSPR: ' num2str(MSE_LUD_GSPR)])
aligned_rots_GSPR = aligned_rots_GSPR_LUD;

scale_factor = iterative_LS_est_scales' * valid_scales / (iterative_LS_est_scales' * iterative_LS_est_scales);
    centered_projections = cryo_addshifts(valid_projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
     [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);
 [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
     [est_shifts,~]=cryo_estimate_shifts(pf,aligned_rots_GSPR,max_shift,shift_step,10000,[],0);
     iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
 centered_projections = cryo_addshifts(valid_projections,-iterative_LS_est_shifts);
     unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));

img_size_GSPR = size(unscaled_centered_projections,1);

    [ v_GSPR, v_b, kernel ,err, iter, flag] = recon3d_firm( unscaled_centered_projections,...
aligned_rots_GSPR,[], 1e-6, 200, zeros(img_size_GSPR,img_size_GSPR,img_size_GSPR));
v_GSPR=real(v_GSPR); v_GSPR(v_GSPR<0)=0;
            save(['simulation_est\flat_vol_estimation_LUD_scale' num2str(max_scale)],...
     'v_GSPR', 'unscaled_centered_projections','iterative_LS_est_shifts','iterative_LS_est_scales',...
     'scale_factor','aligned_rots_GSPR','ref_q_GSPR','est_inv_rots_LUD_GSPR','projections','scales','ref_shifts')
%% Register reconstruction
if register

temp_vol=v_GSPR;
    padlen = (size(temp_vol,1) - size(vol1,1))/2;
    p = padarray(vol1,[padlen,padlen,padlen]);
    corr = fftshift(ifftn(fftn(ifftshift(temp_vol)).*conj(fftn(ifftshift(p))),'symmetric'));
    [~,midx]=max(corr(:));
    [x1,y1,z1] = ind2sub(size(p),midx);
    x2 = size(vol1) - 1;
    eps = 999;
%     for i=linspace(1,1.1,10)
%         temp_vol=imresize3(vol,i);

    for s1=-10:10
        for s2=-10:10
            for s3=-10:10
                xs = round([x1,y1,z1]-x2/2)+[s1,s2,s3];
                xe = round([x1,y1,z1]+x2/2)+[s1,s2,s3];
                if min([xs,xe])<1 || max([xs,xe])>size(temp_vol,1)
                    continue;
                end
                v_new_registered = temp_vol(xs(1):xe(1),xs(2):xe(2),xs(3):xe(3));

                new_eps = (1/sum(vol(:))) * sum(abs(vol1(:)-v_new_registered(:)));
                if new_eps<eps
                    v_GSPR_registered = v_new_registered;
                end   
                eps = min(eps, new_eps);
            end
        end
    end
    
       
   
%     end



            save(['simulation_est\flat_vol_estimation_LUD_scale' num2str(max_scale)],...
      'v_GSPR_registered', 'eps','-append')
end