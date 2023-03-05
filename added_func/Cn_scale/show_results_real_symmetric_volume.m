%% show results
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd('../../');
initpath;
cd(myPath);
%% load results
load('results/sym_flat_real_vol_est.mat')
final_vol(final_vol<max(final_vol(:))*0.01)=0;
final_v = zeros(275,275,275);
final_v(:,:,276/2-18:276/2+17) = final_vol;
%figure;vol3d('cdata',astra_vol);
[MSE,PSNR,recon_proj] = MSEvol(final_v,unscaled_centered_projections,rots);
recon_proj(recon_proj<0)=0;
mean(PSNR)
%figure;compareproj(unscaled_centered_projections,recon_proj,1,4,{},MSE,0);
figure;viewstack(unscaled_centered_projections,4,4);
figure;viewstack(recon_proj,4,4);
% 
% for i=1:5
%     figure
%     imagesc(unscaled_centered_projections(:,:,i));
%         axis off
% axis image
% figure
%     imagesc(recon_proj(:,:,i));
%     axis off
% axis image
% end
%%