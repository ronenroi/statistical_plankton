%% new method
clear;
close all
initpath;
%%
p=parpool(8);
vmax=linspace(0,0.15,6);
n_theta=180;
shift_step = 1;
max_shift = ceil(2*sqrt(2)*15) ;
clusteringThresh = 1;


for i=1:6
load(['deform_vol_projection_' num2str(vmax(i)) '.mat'])
projections = data.projections;
mask_radius = round(0.5*size(projections,1));
n_r = round(0.5*size(projections,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [ unscaled_centered_projections, iterative_LS_est_scales, iter ,...
    clstack, corrstack, iterative_LS_est_shifts, ref_q,ref_shifts,ref_scales] = estimate_scales_shifts(projections, max_shift,...
                                                      max_scale, mask_radius, n_r, n_theta, clusteringThresh, data.ref_q,data.ref_shifts,data.ref_scales);
[est_inv_rots_LUD_GSPR] = est_orientations_LUD(clstack_GSPR, n_theta); % Least Unsquared Deviations
[MSE_LUD_GSPR(kk), ~, ~, aligned_rots_GSPR_LUD] = check_MSE(est_inv_rots_LUD_GSPR, ref_q_GSPR);
disp(['MSE LUD GSPR: ' num2str(MSE_LUD_GSPR(kk))])
aligned_rots_GSPR = aligned_rots_GSPR_LUD;


%%
scale_factor = iterative_LS_est_scales' * valid_scales / (iterative_LS_est_scales' * iterative_LS_est_scales);
    centered_projections = cryo_addshifts(valid_projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
     [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);
 [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
     [est_shifts,~]=cryo_estimate_shifts(pf,aligned_rots_GSPR,max_shift,shift_step,10000,[],0);
     iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
 centered_projections = cryo_addshifts(valid_projections,-iterative_LS_est_shifts);
     unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
save(['3DPOP_deform_results'  num2str(vmax(i)) '.mat'])
   
end