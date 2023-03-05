function cn_symmetry_vol_est(recon_mat_fname,projections,n_symm,...
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

save (recon_mat_fname,'projections','n_symm',...
    'n_theta','n_r','max_shift','shift_step','max_scale','mask_radius','inplane_rot_res','n_Points_sphere','-append')

initstate; 
log_message('symmetry class is C%d',n_symm);

% figure; viewstack(projs,5,5);
if mask_radius > 0
    log_message('Masking projections. Masking radius is %d pixels',mask_radius);
    masked_projs = mask_fuzzy(projections, mask_radius);
else
    masked_projs = projections;
end
log_message('Computing the polar Fourier transform of projections');
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,'double');
log_message('gassian filtering images');
npf = gaussian_filter_imgs(npf);
save(recon_mat_fname,'npf','-append');

log_message('determining third rows outer product using maximum likelihood');
log_message('Maximum shift is %d pixels',max_shift);
log_message('Shift step is %d pixels',shift_step);
log_message('Maximum scale is %d pixels',max_scale);



% if(n_symm==3 || n_symm==4)
% %     max_shift_1d  = ceil(2*sqrt(2)*max_shift); % TODO: is this needed? if so, where?
%     is_remove_non_rank1 = true;
%     non_rank1_remov_percent = 0.25;
%     [vijs,viis,npf,masked_projs] = compute_third_row_outer_prod_c34(n_symm,npf,max_shift,shift_step,recon_mat_fname,...
%         masked_projs,verbose,is_remove_non_rank1,non_rank1_remov_percent);
% else
%     [vijs,viis] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,cache_file_name,verbose);
% end
%% test

if 1
 % [vijs,viis] = compute_third_row_outer_prod_cn(npf,n_symm,max_shift,shift_step,[],0,gt_q);
[vijs,viis] = compute_third_row_outer_prod_cn_scales_shifts_v2(npf,n_symm,max_shift,shift_step,max_scale,Nscales,inplane_rot_res,n_Points_sphere,gt_q);
log_message('Saving third rows outer prods');
    save(recon_mat_fname,'vijs','viis','-append');
else
        load(recon_mat_fname,'vijs','viis');
end

[vijs,viis] = global_sync_J(vijs,viis);
    log_message('Saving third rows outer prods after global_sync');
    save(recon_mat_fname,'vijs','viis','-append');

vis  = estimate_third_rows(vijs,viis,n_symm);
    log_message('Saving third rows under');
    save(recon_mat_fname,'vis','-append');
%         if max_scale == 1
%rots_noscale = estimate_inplane_rotations(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step);
%[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots_noscale,n_symm,n_theta,gt_q);

   
%         else
            rots = estimate_inplane_rotations_shift_scale(npf,vis,n_symm,inplane_rot_res,max_shift,shift_step,max_scale,Nscales);
[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,gt_q);

%         end

        
            log_message('Saving estimated rotation matrices');
    save(recon_mat_fname,'rots','-append');


log_message('Reconstructing abinitio volume');
[sym_estimated_Vol,est_rotations,est_shifts,est_scales] = reconstruct_cn(masked_projs,rots,n_symm,n_r,n_theta,max_shift,shift_step,max_scale,Nscales);
save(recon_mat_fname,'sym_estimated_Vol','est_rotations','est_shifts','-append');
save(recon_mat_fname,'est_scales','-append');
[estR,estdx,est_vol_register,reflect ] = registeration( gt_vol,sym_estimated_Vol );
log_message('Saving estimated vol, rotations and shifts');
save_vols(est_vol_register,recon_mat_fname,n_symm);



end
