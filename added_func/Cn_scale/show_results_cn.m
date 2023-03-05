close all
myPath = fileparts(mfilename('fullpath'));
cd (myPath);
cd('../../');
initpath;
cd('../');
astradir = cd;
addpath(genpath(fullfile([astradir,'/astra-1.9.0.dev11'])));
cd(myPath);
%%

vol_vec{1} = gt_vol;
vol_vec{1}(gt_vol<0)=0;
% vol_vec{2} = recon_gt;
% vol_vec{2}(recon_gt<0) = 0;
 vol_vec{2} = vol_filt_gn;
 vol_vec{2}(vol_filt_gn<0) = 0;
% vol_vec{4} = vol_no_filt;
% vol_vec{4}(vol_no_filt<0) = 0;
% [estR,estdx,est_vol_register,reflect ] = registeration( gt_vol,sym_estimated_Vol );
% [estR,estdx,est_vol_register,reflect]=cryo_align_densities(gt_vol,est_vol_register);
% vol_vec{2} = est_vol_register;
% vol_vec{2}(est_vol_register<0) = 0;
% 
% [estR,estdx,est_vol_register,reflect ] = registeration( gt_vol,vol_filt );
% [estR,estdx,est_vol_register,reflect]=cryo_align_densities(gt_vol,est_vol_register);
% vol_vec{3} = est_vol_register;
% vol_vec{3}(est_vol_register<0) = 0;
% 
% [estR,estdx,est_vol_register,reflect ] = registeration( gt_vol,vol_no_filt );
% [estR,estdx,est_vol_register,reflect]=cryo_align_densities(gt_vol,est_vol_register);
% vol_vec{4} = est_vol_register;
% vol_vec{4}(est_vol_register<0) = 0;

[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,gt_q(:,valid_eq));
% save('sym_vol_est_shifts_scales1.mat','vol_vec','-append')
titles_vec{1} = 'Ground truth';
titles_vec{2} = 'Reconstruction volume 1';
titles_vec{3} = 'Reconstruction volume 2';
titles_vec{4} = 'Reconstruction volume 3';

err=1;

for i=1:length(vol_vec)
    figure;
    curr_vol = vol_vec{i};
    
    thr = max(curr_vol(:)) * 30 / 100;
    err = norm(curr_vol(:)-vol_vec{1}(:))/norm(vol_vec{1}(:))
    curr_vol(curr_vol < thr) = 0;
    %subplot(2,ceil(length(vol_vec)/2),i);
    vol3d('cdata',curr_vol);
    %title({titles_vec{i},['relative error = ' num2str(err)]});
%    title({titles_vec{i}});

     colormap(hsv);
    %caxis([0 max(vol_vec{i}(:))])
    % alphamap(alpha);
    axis equal off
    set(gcf, 'color', 'w');
    view(3);
end

%% projections comparison
n=5;                              % projection number
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
    test_projections{i} = permute(proj,[2 1 3]);   % transpose each image
    lim_max = max(max(test_projections{i}(:)),lim_max);
    lim_min = min(min(test_projections{i}(:)),lim_min);
    
end
% show projections
figure
colormap(gray)
for i=1:length(test_projections)
    curr_proj = test_projections{i};
    for j=1:size(curr_proj,3)
        subplot(length(test_projections),size(curr_proj,3),(i-1)*size(curr_proj,3)+j)
        curr_im = curr_proj(:,:,j);
        imagesc(curr_im);
        axis off
        axis image
        title({titles_vec{i},['proj. #' int2str(j)]});
    end
    caxis([lim_min lim_max])
    
end
