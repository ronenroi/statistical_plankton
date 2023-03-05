%% show results
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd('../../');
initpath;
cd(myPath);
%% load results
load('est. results/jules_est/daphnia_est_registered_new_met.mat')
%  
%figure;vol3d('cdata',astra_vol);

    img = unscaled_centered_projections;
    rots = rotations;
    [vol_astra] = astra_reconstruction(img,rots,[165 165 165],1000);
    vol_astra(vol_astra<0.3*max(vol_astra(:)))=0;
[MSE_SEM(i),PSNR_SEM(i), recon_proj] = MSEvol(vol_astra,unscaled_centered_projections,rotations);    


%%
load('est. results/jules_est/daphnia_est_registered_new_met.mat')
%  
%figure;vol3d('cdata',astra_vol);

for i=1:size(unscaled_centered_projections,3)
    idx = 1:size(unscaled_centered_projections,3);
    idx(i)=[];
    img = unscaled_centered_projections(:,:,idx);
    rots = rotations(:,:,idx);
    [vol_astra] = astra_reconstruction(img,rots,[165 165 165],1000);
    vol_astra(vol_astra<0.3*max(vol_astra(:)))=0;
[MSE_SEM(i),PSNR_SEM(i), recon_proj] = MSEvol(vol_astra,unscaled_centered_projections(:,:,i),rotations(:,:,i));    
end
save('psnrSEM_daphnia','PSNR_SEM');
%%

load('ICCP\daphnia compare\3DPOP_daphnia_1.mat')
for i=1:size(unscaled_centered_projections,3)
 idx = 1:size(unscaled_centered_projections,3);
    idx(i)=[];
    img = unscaled_centered_projections(:,:,idx);
    rots = est_inv_rots_LUD_GSPR(:,:,idx);
    [vol_astra] = astra_reconstruction(img,rots,[165 165 165],1000);
    vol_astra(vol_astra<0.3*max(vol_astra(:)))=0;
[MSE_pop,PSNR_pop(i), recon_proj_pop] = MSEvol(vol_astra,unscaled_centered_projections,est_inv_rots_LUD_GSPR);
save('psnrPOP_daphnia','PSNR_pop');

end
figure;viewstack(recon_proj_pop,4,4)
figure;viewstack(unscaled_centered_projections,4,4)