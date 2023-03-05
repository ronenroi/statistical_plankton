%% show results
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd('../../');
initpath;
cd(myPath);
%% load results
% load('jules_est/vol_est_star_processed.mat')
%% MSE Aviad
load('est. results/recon_Pyramimonas_longicauda_aviad1.mat')

 v(v<max(v(:))*0.1)=0;
[MSE_aviad,PSNR_aviad,recon_proj_aviad] = MSEvol(v,unscaled_projections,est_inv_rots_LUD);
 mean(PSNR_aviad)
load('est. results/recon_Pyramimonas_longicauda_new_met.mat')
 vol_astra(vol_astra<max(vol_astra(:))*0.1) = 0;
[MSE_roi,PSNR_roi, recon_proj_roi] = MSEvol(vol_astra,unscaled_centered_projections,rotations);
mean(PSNR_roi)
recon_proj_aviad(recon_proj_aviad<0)=0;
recon_proj_roi(recon_proj_roi<0)=0;

figure;compareproj(unscaled_projections,recon_proj_aviad,2,4,'Upper img - Original Input Projections; Lower img - Reconstructed Volume Projections Aviad',MSE_aviad);
figure;compareproj(unscaled_centered_projections,recon_proj_roi,2,4,'Upper img - Original Input Projections; Lower img - Reconstructed Volume Projections New code',MSE_roi);


%% 3D plots


vol_vec{1} = v;
vol_vec{1}(v<0)=0;
if ~exist('vol_aligned','var')
[resA,vol_aligned,h]=cryo_compare_volumes(est_vol_new,v);
save('est. results/recon_Pyramimonas_longicauda_new_met.mat','vol_aligned','-append')
end
vol_vec{2} = vol_aligned;
vol_vec{2}(vol_aligned<0)=0;


titles_vec{1} = 'Reconstruction Pyramimonas longicauda (Avaid)';
titles_vec{2} = 'Reconstruction Pyramimonas longicauda (Roi)';

figure;
for i=1:length(vol_vec)
    figure;

    curr_vol = vol_vec{i};
    thr = max(vol_vec{i}(:)) * 30 / 100;
    curr_vol(curr_vol < thr) = 0;
   % subplot(2,1,i);
    vol3d('cdata',curr_vol);
     title(titles_vec{i});
%     colormap(gray(256));
    colormap(jet(256));
    caxis([0 max(vol_vec{1}(:))])
      alphamap([0:0.01:0.3])
     axis  normal off
    set(gcf, 'color', 'w');
    view(3);  
end




