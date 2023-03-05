%% show results
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd('../../');
initpath;
cd(myPath);
%% load results
load('simulation_est/simulation_est_GSPR_registered_LUD_scale.mat')
% load('simulation_est/simulation_estimation_registered_LUD_scale1.mat')
load('simulation_est/simulation_est_registered_new_met_final.mat')
load ('simulation_est/gt_recon.mat')
load volume
%% 3D plots
% alpha = exp(linspace(-4,4,255));
% alpha = alpha - min(alpha(:));
% alpha = alpha / max(alpha(:));
% vol{1} = data.ref_vol;
vol_vec{1} = vol;
vol_vec{1}(vol<0)=0;
vol_vec{2} = est_vol_gt_register;
vol_vec{2}(est_vol_gt_register<0) = 0;
vol_vec{3} = v_GSPR_registered;
vol_vec{3}(vol_vec{3}<0)=0;
vol_vec{4} = est_vol_new_register;
vol_vec{4}(vol_vec{4}<0)=0;

titles_vec{1} = 'Ground truth';
titles_vec{2} = {'Reconstruction with GT',' scale-shift-rotations'};
titles_vec{3} = '3D-POP';
titles_vec{4} = '3D-PMPO (ours)';


for i=1:length(vol_vec)
    figure;
    curr_vol = vol_vec{i};
    thr = max(vol_vec{1}(:)) * 20 / 100;
    err(i) = norm(curr_vol(:)-vol_vec{1}(:))/norm(vol_vec{1}(:));
    curr_vol(curr_vol < thr) = 0;
%     subplot(1,length(vol_vec),i);
    vol3d('cdata',curr_vol);
    title(titles_vec{i});
    colormap(jet(256));
    caxis([0 max(vol_vec{1}(:))])
    % alphamap(alpha);
    axis  normal off
    set(gcf, 'color', 'w');
    view(3);
end

%% projections comparison
n=4;                              % projection number
q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
lim_max = -inf;
lim_min = inf;

for k = 1:n
    quat = q(:, k);
    q(:, k) = quat / norm(quat);
end
for i=1:length(vol_vec)
    curr_vol = vol_vec{i};
    proj = cryo_project(curr_vol,q_to_rot(q),size(curr_vol,1),'single'); % generate projecitons
    proj(proj<0)=0;
    test_projections{i} = permute(proj,[2 1 3]);   % transpose each image
    lim_max = max(max(test_projections{i}(:)),lim_max);
    lim_min = min(min(test_projections{i}(:)),lim_min);

end
%% show projections

    colormap(jet(256));
for i=1:length(test_projections)
    figure
    curr_proj = test_projections{i};
    for j=1:size(curr_proj,3)
        subplot(1,size(curr_proj,3),j)
        curr_im = curr_proj(:,:,j);
        imagesc(curr_im);
        axis off
        axis image
       % title(titles_vec{i});
    end
        caxis([lim_min lim_max])

end


    
