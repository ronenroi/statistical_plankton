function[gt_vol,projections,inv_rot_matrices,q_ref,ref_shifts,ref_scales,n_symm] = generate_symetric_projections(n_Images,max_scale)

n_symm = 8;
%% create n_symm volume
im_size = 129;
gt_vol = zeros(im_size,im_size,im_size);
x0 = round(im_size/2);
[X,Y] = meshgrid(linspace(-x0,x0,im_size),linspace(x0,-x0,im_size));
zmax = 5;
for z=0:zmax
plane=zeros(im_size,im_size);
m = double((X.^2+Y.^2).^0.5 < x0 * 0.8 * cos(z/zmax*pi/2)* sin(atan(Y./X)*n_symm));
 m = m .* (1+z)/zmax;
 cir =double(((X-x0 * 0.3*sin(atan(Y./X)*n_symm) ).^2+(Y-x0 * 0.3*sin(atan(Y./X)*n_symm) ).^2).^0.5 < 5 *sin(atan(Y./X)*n_symm) );
 m = m + cir * (zmax - z)*2;
for i=0:n_symm-1
plane = plane + imrotate(m,45*i,'bilinear','crop');
end
gt_vol(:,:,z+x0) = plane;
if z >0
    gt_vol(:,:,-z+x0) = plane;
end
end
figure;vol3d('Cdata',gt_vol);view(3);
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
        projections=cryo_project(gt_vol,q_to_rot (q_ref),size(gt_vol,1),'single'); % generate phantom projecitons

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

    % Add shifts to projections
    max_shift = 10;
    [projections, ref_shifts] = cryo_addshifts(scaled_projections,[],max_shift,1);
    %projections=scaled_projections; ref_shifts=0;
    projections(projections<0)=0;
    max_proj = max(projections(:));
    max_well = 20000;
    dark_current = 10 * 1e-6 * ones(size(projections(:,:,1)));
    % Add noise 
    for i = 1:n_Images
        projections(:,:,i) = max_proj/(max_well*0.5*1e-6) * ...
            (imnoise(single(dark_current),'poisson') + imnoise(imnoise(0.5*max_well*1e-6*single(projections(:,:,i)/max_proj), 'poisson'),'speckle', 1e-6));
    end

%figure;viewstack(projections,5,5);   % Display the proejctions.

end

