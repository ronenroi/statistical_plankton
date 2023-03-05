function cn_symmetry_flat_vol_est(recon_mat_fname,projections,n_symm,...
    n_theta,n_r,max_shift,shift_step,max_scale,Nscales,mask_radius,inplane_rot_res,n_Points_sphere,gt_q,gt_vol)
if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_r','var')
    n_r = round(0.75*size(projections,1));
end

if ~exist('max_shift','var')
    max_shift = 15;
end
if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('mask_radius','var')
    mask_radius = round(0.75*size(projections,1));
    
end

if ~exist('inplane_rot_res','var')
    inplane_rot_res = 1;
end
if ~exist('n_Points_sphere','var')
    n_Points_sphere = 100;
end
if ~exist('max_scale','var')
    max_scale = 1;
end

%save (recon_mat_fname,'projections','n_symm','n_theta','n_r','max_shift','shift_step','max_scale','mask_radius','inplane_rot_res','n_Points_sphere');

initstate;
log_message('symmetry class is C%d',n_symm);

% figure; viewstack(projs,5,5);
%% Orientation assigment
thr = 0.05;

[projectionDir,ellipseMeanstack,iterative_LS_est_scales] = ellipse_estimation(projections,thr);
iterative_LS_est_shifts=[ellipseMeanstack(:,2) ellipseMeanstack(:,1)];
   centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./iterative_LS_est_scales);
   



%%
clusteringThresh=0.15;
discard_projections = [];
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
%%
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(unscaled_centered_projections, mask_radius);
    log_message('Computing the polar Fourier transform of projections');
    [npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');
    log_message('gassian filtering images');
    npf_filt = gaussian_filter_imgs(npf);
    %save(recon_mat_fname,'npf','-append');
    %%
    
    %%
    log_message('determining third rows outer product using maximum likelihood');
    log_message('Maximum shift is %d pixels',max_shift);
    log_message('Shift step is %d pixels',shift_step);
    log_message('Maximum scale is %d pixels',max_scale);
    [vijs,viis] = compute_third_row_outer_prod_cn_scales_shifts_flat(npf_filt,n_symm,max_shift,shift_step,max_scale,Nscales,inplane_rot_res,n_Points_sphere,projectionDir);
    
    %log_message('Saving third rows outer prods');
    save(recon_mat_fname,'vijs','viis','-append');
    [vijs,viis] = global_sync_J(vijs,viis);
    vis  = estimate_third_rows(vijs,viis,n_symm);
    [rots,~,corrstack] = estimate_inplane_rotations_shift_scale(npf_filt,vis,n_symm,inplane_rot_res,max_shift,shift_step,max_scale,Nscales);
    save(recon_mat_fname,'-append');
    COV = corrstack + corrstack' + eye(size(corrstack,1));
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1)<clusteringThresh;
    if any(not(valid_eq))
        discard_projections = cat(3,discard_projections,projections(:,:,not(valid_eq)));
        projections = projections(:,:,valid_eq);
        
%         iterative_LS_est_scales = iterative_LS_est_scales(valid_eq);
%         iterative_LS_est_shifts = iterative_LS_est_shifts(valid_eq,:);
        npf_filt=npf_filt(:,:,valid_eq);
        npf = npf(:,:,valid_eq);
        rots = rots(:,:,valid_eq);
        iterative_LS_est_shifts=iterative_LS_est_shifts(valid_eq,:);
        iterative_LS_est_scales=iterative_LS_est_scales(valid_eq);
        if exist('gt_q') 
            gt_q = gt_q(:,valid_eq);
        end
        fprintf('Removed %i projections \n',sum(not(valid_eq)));
    end
    
    
    
    [est_shifts,est_scales] = max_estimate_shifts_scales_cn(npf,rots,n_symm,max_shift,shift_step,max_scale,Nscales);
    iterative_LS_est_scales = iterative_LS_est_scales.*est_scales;
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./iterative_LS_est_scales);
    
    
    
    max_scale = max(abs(1-est_scales))+1;
    max_shift = max(sqrt(sum(abs(est_shifts).^2,2))) / 2;
    save(recon_mat_fname,'-append');
log_message('Reconstructing volume');
[sym_estimated_Vol_gn,est_rotations,est_shifts_cn,est_scales_cn] = reconstruct_cn(unscaled_centered_projections,rots,n_symm,n_r,n_theta,max_shift,shift_step,max_scale,Nscales);
% [sym_estimated_Vol, ~, ~ ,~, ~, ~] = recon3d_firm_parallel( unscaled_centered_projections,rots,[], 1e-6, 200, zeros(size(unscaled_centered_projections,1),size(unscaled_centered_projections,1),size(unscaled_centered_projections,1)));
save(recon_mat_fname,'-append');
if exist('gt_vol','var')
[~,~,sym_estimated_Vol_gn_reg,~ ] = registeration( gt_vol,sym_estimated_Vol_gn );
[~,~,sym_estimated_Vol_reg,~ ] = registeration( gt_vol,sym_estimated_Vol );
save(recon_mat_fname,'-append');
end
log_message('Saving all data');
vol_filt_gn    = make_vol_cn(sym_estimated_Vol_gn,n_symm,true);
vol_no_filt_gn = make_vol_cn(sym_estimated_Vol_gn,n_symm,false);
% vol_filt    = make_vol_cn(sym_estimated_Vol_reg,n_symm,true);
% vol_no_filt = make_vol_cn(sym_estimated_Vol_reg,n_symm,false);
save(recon_mat_fname,'-append');

end
