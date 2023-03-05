%% new method
clear;
close all
initpath;
for i=1:6
    vmax=linspace(0,0.15,6);
    load(['SEM_deform_results'  num2str(vmax(i)) '.mat'],'iterative_LS_est_scales','iterative_LS_est_shifts','data','aligned_rots','unscaled_centered_projections');
      img_size = size(unscaled_centered_projections,1);
 [ v1, ~, ~ ,~, ~, ~] = recon3d_firm( unscaled_centered_projections,...
     aligned_rots,[], 1e-6, 200, zeros(img_size,img_size,img_size));
 est_vol_new=real(v1);
 est_vol_new(est_vol_new<0)=0;
     save(['SEM_deform_vols' num2str(vmax(i)) '.mat']);
     [estR,estdx,est_vol_register,reflect,estscale,eps_beta_SEM ] = registeration( data.ref_vol,est_vol_new );
        save(['SEM_deform_vols' num2str(vmax(i)) '.mat']);
     clear all;
end
