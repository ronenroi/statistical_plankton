function[gt_vol,projections,inv_rot_matrices,q_ref,ref_shifts,ref_scales] = generate_circle_projections(n_Images,max_shift,max_scale,AddNoise)
%% create circle volume
im_size = 129;
gt_vol = zeros(im_size,im_size,im_size);
x0 = round(im_size/2);
[X,Y] = meshgrid(linspace(-x0,x0,im_size),linspace(x0,-x0,im_size));
plane = double((X.^2+Y.^2).^0.5 < x0 * 0.8 );
gt_vol(:,:,x0) = plane;
%figure;vol3d('Cdata',gt_vol);view(3); axis('equal')
%% create projections

    inv_rot_matrices = zeros(3, 3, n_Images);
    q_ref = qrand(n_Images);                       % generate rotations as quaternions
    % find inverse rotation matrices
    for k = 1:n_Images
        quat = q_ref(:, k);
        q_ref(:, k) = quat/norm(quat);
        rot = q_to_rot(q_ref(:, k));
        inv_rot_matrices(:, :, k) = rot';
    end
%     gt_theta = rand(1,n_Images)*90;
%     gt_phi = rand(1,n_Images)*90;
%    for k = 1:n_Images
%         rot = euler_to_rot([0;gt_theta(k);gt_phi(k)]);
%        rot_matrices(:, :, k) = rot;
%     end
        projections=cryo_project(gt_vol,q_to_rot(q_ref),size(gt_vol,1),'single'); % generate phantom projecitons

    % The following step is ESSENTIAL before calling FIRM reconstruction on the
    % projections.
    projections=permute(projections,[2 1 3]);   % transpose each image
    
    % Pad projections 
    pad_factor = 1.15*max_scale;
    
    imsize = size(projections,1);
    padded_size = floor(pad_factor*[imsize, imsize]);

    padded_projections = padarray(projections, ceil((padded_size-imsize)/2), 'pre');
    padded_projections = padarray(padded_projections, ceil((padded_size-imsize)/2), 'post');

    ref_scales = (max_scale - 1/max_scale)*rand(n_Images,1) + 1/max_scale;
    scaled_projections = scale_projections(padded_projections, ref_scales);
    [projections, ref_shifts] = cryo_addshifts(scaled_projections,[],max_shift,1);
    projections(projections<0)=0;
    if AddNoise
    % Add shifts to projections

    max_proj = max(projections(:));
    max_well = 20000;
    dark_current = 10 * 1e-6 * ones(size(projections(:,:,1)));
    % Add noise 
    for i = 1:n_Images
        projections(:,:,i) = max_proj/(max_well*0.5*1e-6) * ...
            (imnoise(single(dark_current),'poisson') + imnoise(imnoise(0.5*max_well*1e-6*single(projections(:,:,i)/max_proj), 'poisson'),'speckle', 1e-6));
    end
    end
figure;viewstack(projections,5,5);   % Display the proejctions.

end


