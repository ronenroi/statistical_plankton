%%
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd (myPath);
cd ('../../../');
initpath;
cd(myPath);
%%
n_Images = 50;
max_scale = 1.5;%1 for no scaling
max_shift = 10; %0 for no shift
AddNoise = true;
 dataFolder = '../../../../data sets/jules images/star/processed/usable'; projectionType = 'DarkField';  
 projections = func_loadPlanktonData(dataFolder, projectionType);
n_Images = size(projections,3);

[gt_vol,projections,inv_rot_matrices,gt_q,gt_shifts,gt_scales] = generate_symetric_projections(n_Images,max_scale);
%[gt_vol,projections,gt_theta,gt_q,gt_shifts,gt_scales] = generate_circle_projections(n_Images,max_shift,max_scale,AddNoise);
%% your code
thr = 0.05;
[ellipsestack,centers,propsstack] = imfindellipse(projections,thr,3,1 );
    MinorAxisLength=[propsstack.MinorAxisLength];
    MajorAxisLength=[propsstack.MajorAxisLength];
for i=1:n_Images
    
    ellipse = ellipsestack(:,:,i);
    im_mask = bwmorph(ellipse, 'bridge');
    im_mask = imfill(im_mask,'holes');
    im_mask_stack(:,:,i)=im_mask;
    [nx, ny, nz] = size(ellipse);
    [X,Y] = meshgrid(linspace(1,nx,nx),linspace(1,ny,ny));
    ellipseSize = sum(im_mask(:));
    ellipseSize1 = sum(sum(ellipsestack(:,:,i)));
    ellipseMean = 1/ellipseSize * [sum(sum(X.*im_mask)),sum(sum(Y.*im_mask))];
    Covxx =  1/ellipseSize * sum(sum(((X - ellipseMean(1)).^2 ).*im_mask));
    Covxx1 = 1/ellipseSize1 * sum(sum(((X - centers(i,1)).^2 ).*ellipsestack(:,:,i)));
    Covyy =  1/ellipseSize * sum(sum(((Y - ellipseMean(2)).^2 ).*im_mask));
    Covyy1 = 1/ellipseSize1 * sum(sum(((Y - centers(i,2)).^2 ).*ellipsestack(:,:,i)));
    Covxy =  1/ellipseSize * sum(sum(((X - ellipseMean(1)).*(Y - ellipseMean(2))).*im_mask));
    ellipseCov = [Covxx, Covxy; Covxy, Covyy];
    [V,D] = eig(ellipseCov);
    a=2 * sqrt(diag(D));
    a1=[MinorAxisLength(i)/2, MajorAxisLength(i)/2];
    ellipseMeanstack(i,:) = ellipseMean;
    projectionDir(:,i) = 1/a(2) * [sqrt(a(2)^2 -4 * Covxx), sqrt(a(2)^2 -4 * Covyy), a(1)];
    projectionDir1(:,i) = 1/a1(2) * [sqrt(a1(2)^2 -4 * Covxx1), sqrt(a1(2)^2 -4 * Covyy1), a1(1)];

    if sign(projectionDir(1,i)*projectionDir(2,i)) ~= -sign(Covxy)
        projectionDir(1,i)=-projectionDir(1,i);
    end
    est_scales(i) = a(2);
end
%% test center and axis
shift1 = mean(mean((centers-112 - [gt_shifts(:,2),gt_shifts(:,1)]).^2))
shift2 = mean(mean((ellipseMeanstack-112 - [gt_shifts(:,2),gt_shifts(:,1)]).^2))
scale1 = mean((([propsstack.MajorAxisLength]/min([propsstack.MajorAxisLength]))' - gt_scales/min(gt_scales(:))).^2)
scale2 = mean(((est_scales/min(est_scales))' - gt_scales/min(gt_scales(:))).^2)
rot1 = mean(mean((abs(projectionDir1) - abs(squeeze(inv_rot_matrices(3,:,:)))).^2))
rot2 = mean(mean((abs(projectionDir) - abs(squeeze(inv_rot_matrices(3,:,:)))).^2))
%%
shift2 = mean(mean((ellipseMeanstack-112 - [gt_shifts(:,2),gt_shifts(:,1)]).^2))
alpha = est_scales * gt_scales / (est_scales' * est_scales);
scale2 = mean(((est_scales/(alpha))' - gt_scales).^2)
rot2 = mean(mean((abs(projectionDir) - abs(squeeze(inv_rot_matrices(3,:,:)))).^2))


%% 4 option for 3rd row of inv rot +a1 +a2 +a3 / +a1 +a2 -a3 / -a1 -a2 +a3 / -a1 -a2 -a3