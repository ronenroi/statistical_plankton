%% new method
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../');
initpath;
cd(myPath);
%%
load('deform_vol_projection_0.06_4.mat', 'data')
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = 15;


%%
n_theta=180;
shift_step = 1;
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));
SEM_worker(projections,'SEM_def0.06_4.mat',1,...
    n_theta,n_r,max_shift,shift_step,max_scale,data);
clear
%%
load('deform_vol_projection_0.06_5.mat', 'data')
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = 15;

n_theta=180;
shift_step = 1;
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));
SEM_worker(projections,'SEM_def0.06_5.mat',1,...
    n_theta,n_r,max_shift,shift_step,max_scale,data);
clear
%%

for i=1:5
    
    load(['deform_vol_projection_0.06_' num2str((i)) '.mat'])
    projections = data.projections;
    mask_radius = round(0.75*size(projections,1));
    n_r = round(0.75*size(projections,1));
    n_theta=180;
    shift_step = 1;
    max_shift = ceil(2*sqrt(2)*15) ;
    clusteringThresh = 1;
    max_scale = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ projections, iterative_LS_est_scales, iter ,...
        clstack, corrstack, iterative_LS_est_shifts, ref_q,ref_shifts,ref_scales] = estimate_scales_shifts(projections, max_shift,...
        max_scale, mask_radius, n_r, n_theta, clusteringThresh, data.ref_q,data.ref_shifts,data.ref_scales);
    save(['3DPOP_deform_results_006_'  num2str((i)) '.mat'])
    [est_inv_rots_LUD_GSPR] = est_orientations_LUD(clstack, n_theta); % Least Unsquared Deviations
    [MSE_LUD, ~, ~, aligned_rots_GSPR_LUD] = check_MSE(est_inv_rots_LUD_GSPR, ref_q);
    aligned_rots_GSPR = aligned_rots_GSPR_LUD;
    
    
    %%
    scale_factor = iterative_LS_est_scales' * ref_scales / (iterative_LS_est_scales' * iterative_LS_est_scales);
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
    [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);
    [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
    [est_shifts,~]=cryo_estimate_shifts(pf,aligned_rots_GSPR,max_shift,shift_step,10000,[],0);
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
    %  [vol_astra] = astra_reconstruction(unscaled_centered_projections,aligned_rots_GSPR,[size(projections,1) size(projections,1) size(projections,1)],1000);
    
    save(['3DPOP_deform_results_006_'  num2str((i)) '.mat'])
    
end