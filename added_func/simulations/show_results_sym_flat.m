%% show results
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd('../../');
initpath;
cd(myPath);
%% load results
load('simulation_est/flat_vol_estimation_LUD_scale2.mat')

%% 3D plots
% alpha = exp(linspace(-4,4,255));
% alpha = alpha - min(alpha(:));
% alpha = alpha / max(alpha(:));
% vol{1} = data.ref_vol;
load volume
x0 = round(size(vol,1)/2);
flat_vol = vol(:,:,x0-2:x0+2);

vol_vec{1} = flat_vol;
vol_vec{1}(flat_vol<0)=0;

sym_vol = zeros(size(vol));
sym_vol(x0:end,x0:end,x0:end) = vol(x0:end,x0:end,x0:end) ;
sym_vol = sym_vol + fliplr(sym_vol) ;
sym_vol = sym_vol + flipud(sym_vol) ;
sym_vol = sym_vol + flip(sym_vol,3) ;
vol_vec{2} = sym_vol;
vol_vec{2}(sym_vol<0)=0;

vol_vec{3} = v_GSPR_registered;
vol_vec{3}(v_GSPR_registered<0)=0;
load('simulation_est/sym_vol_estimation_LUD_scale2.mat')

vol_vec{4} = v_GSPR_registered;
vol_vec{3}(v_GSPR_registered<0) = 0;

titles_vec{1} = 'GT flat vol';
titles_vec{2} = 'GT sym vol';
titles_vec{3} = 'flat vol est';
titles_vec{4} = 'sym vol est';



figure;
for i=1:length(vol_vec)
    curr_vol = vol_vec{i};
    thr = max(vol_vec{i}(:)) * 10 / 100;
%     err = norm(curr_vol(:)-vol_vec{1}(:))/norm(vol_vec{1}(:));
    curr_vol(curr_vol < thr) = 0;
    subplot(2,length(vol_vec)/2,i);
    vol3d('cdata',curr_vol);
%     title({titles_vec{i},['relative error = ' num2str(err)]});
     title(titles_vec{i});

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
%% show projections
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


    
