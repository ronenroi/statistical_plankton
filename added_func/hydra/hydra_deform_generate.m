%%create deform hydra
%% new method
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
%% load volume

for kkk = 1:6
load('data_6.mat','ref_vols')
    
vol = ref_vols{kkk};
[nx, ny, nz] = size(vol);
n=100;                              

vmax = 20;
inv_rot_matrices = zeros(3, 3, n);
q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
for k = 1:n
	quat = q(:, k);
	q(:, k) = quat/norm(quat);
	rot = q_to_rot(q(:, k));
	inv_rot_matrices(:, :, k) = rot';
end
 for   j=1:n
[ deformed_vol1 ] = deform(deformed_vol, vmax);

    projections(:,:,j)=cryo_project(deformed_vol,q_to_rot(q(:,j)),size(deformed_vol,1),'single'); % generate phantom projecitons
end
 projections=permute(projections,[2 1 3]);   % transpose each image
max_scale = 2.2;
max_log_scale = log(max_scale);

log_scales = 2*max_log_scale*rand(n,1) - max_log_scale;
log_scales = log_scales - mean(log_scales);
scales = exp(log_scales);
    % pad projections 
    pad_factor = 2;
    
    imsize = size(projections,1);
    padded_size = floor(pad_factor*[imsize, imsize]);

    padded_projections = padarray(projections, ceil((padded_size-imsize)/2), 'pre');
    padded_projections = padarray(padded_projections, ceil((padded_size-imsize)/2), 'post');

    scaled_projections = scale_projections(padded_projections, scales);

    % add shifts to projections
    max_shift = 10;
    [projections, ref_shifts] = cryo_addshifts(scaled_projections,[],max_shift,1);
    %projections=scaled_projections; ref_shifts=0;
    projections(projections<0)=0;
    max_proj = max(projections(:));
    max_well = 20000;
    dark_current = 10 * 1e-6 * ones(size(projections(:,:,1)));
    % add noise 
    for j = 1:n
        projections(:,:,j) = max_proj/(max_well*0.5*1e-6) * ...
            (imnoise(single(dark_current),'poisson') + imnoise(imnoise(0.5*max_well*1e-6*single(projections(:,:,j)/max_proj), 'poisson'),'speckle', 1e-6));
    end

    data.ref_q = q;
    data.ref_scales = scales;
    data.ref_shifts = ref_shifts;
    data.inv_rot_matrices = inv_rot_matrices;
    data.ref_vol = vol;
    data.projections = projections;
    data.max_shift = max_shift;
    data.max_scale = max_scale;

    
	save(['deform_vol_projection1_10_' num2str(kkk) '.mat'],'data');
 clear;


end