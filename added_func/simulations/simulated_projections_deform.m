%% load volume
load volume.mat

max_scale = 2;
[nx, ny, nz] = size(vol);
n=100;                              
inv_rot_matrices = zeros(3, 3, n);
q = qrand(n);                       % generate rotations as quaternions
% find inverse rotation matrices
for k = 1:n
	quat = q(:, k);
	q(:, k) = quat/norm(quat);
	rot = q_to_rot(q(:, k));
	inv_rot_matrices(:, :, k) = rot';
end
scales = (max_scale - 1/max_scale)*rand(n,1) + 1/max_scale;                   
vmax=linspace(0,0.15,6);
for i=1:length(vmax)
 for   j=1:n
[ deformed_vol ] = deform(vol, vmax(i));

    projections(:,:,j)=cryo_project(deformed_vol,q_to_rot(q(:,j)),size(deformed_vol,1),'single'); % generate phantom projecitons
end
    % the following step is essential before calling firm reconstruction on the
    % projections.
    projections=permute(projections,[2 1 3]);   % transpose each image

    % pad projections 
    pad_factor = 1.15*max_scale;
    
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
    data.ref_shifts = shifts;
    data.inv_rot_matrices = inv_rot_matrices;
    data.ref_vol = vol;
    data.projections = projections;
    data.max_shift = max_shift;
    data.max_scale = max_scale;
    
	save(['deform_vol_projection_' num2str(vmax(i))],'data');
end
