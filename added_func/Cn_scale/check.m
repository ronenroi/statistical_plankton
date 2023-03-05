[rot_alligned,err_in_degrees,mse] = analyze_results_cn(rots,n_symm,n_theta,gt_q);

% 
% 
% 
% nImages = size(rots,3);
% g = [cosd(360/n_symm) -sind(360/n_symm) 0;
% sind(360/n_symm)  cosd(360/n_symm) 0;
% 0 				 0  1]; % a rotation of 90 degrees about the z-axis
% for k=1:nImages
% rot = rots(:,:,k);
% for s=0:n_symm-1
% RsCn(:,:,s*nImages+k) = g^s*rot;
% end
% end
% projsCn = [];
% for i=1:n_symm
% projsCn = cat(3,projsCn,projections);
% %     projsCn(:,:,i) = projs;
% end
%%
centered_projections = cryo_addshifts(projections,-est_shifts);
centered_scaled_projections = scale_projections(centered_projections, 1.0./est_scales);
%%
% projsCn=permute(projsCn,[1 3 2]);   % col x n_proj x row
% 
% vectors = rot2astraVec(RsCn);
% vol_geom = astra_create_vol_geom([297,297,297]);
% proj_geom = astra_create_proj_geom('parallel3d_vec', size(projsCn,3), size(projsCn,1), vectors);
% 
% [vol_id, recon_vol] = astra_create_backprojection3d_cuda(projsCn, proj_geom, vol_geom);
% astra_mex_data3d('delete', vol_id);
%%
[ recon_vol, ~, ~ ,~, ~, ~] = recon3d_firm( centered_scaled_projections,rot_alligned,-[], 1e-6, 100, zeros(297,297,297));

%% 3D plots

vol_vec{1} = gt_vol;
vol_vec{1}(gt_vol<0)=0;
vol_vec{2} = recon_vol;
vol_vec{2}(recon_vol<0) = 0;


titles_vec{1} = 'Ground truth';
titles_vec{2} = 'Reconstruction with gt scale-shift';
err=1;

figure;
for i=1:length(vol_vec)
    curr_vol = vol_vec{i};
    thr = max(vol_vec{i}(:)) * 60 / 100;
%     err = norm(curr_vol(:)-vol_vec{1}(:))/norm(vol_vec{1}(:));
    curr_vol(curr_vol < thr) = 0;
    subplot(2,length(vol_vec)/2,i);
    vol3d('cdata',curr_vol);
    title({titles_vec{i},['relative error = ' num2str(err)]});
    colormap(gray(256));
    caxis([0 max(vol_vec{i}(:))])
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
