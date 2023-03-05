function [projectionDir,ellipseMeanstack,est_scales] = ellipse_estimation(projections,thr)
n_Images = size(projections,3);
[ellipsestack,~,~] = imfindellipse(projections,thr,3,1 );
ellipseMeanstack = zeros(n_Images,2);
projectionDir = zeros(3,n_Images);
est_scales = zeros(n_Images,1);
for i=1:n_Images
    
    ellipse = ellipsestack(:,:,i);
    im_mask = bwmorph(ellipse, 'bridge');
    im_mask = imfill(im_mask,'holes');
    %im_mask_stack(:,:,i)=im_mask;
    [nx, ny] = size(ellipse);
    [X,Y] = meshgrid(linspace(1,nx,nx),linspace(1,ny,ny));
    ellipseSize = sum(im_mask(:));
    ellipseMean = 1/ellipseSize * [sum(sum(X.*im_mask)),sum(sum(Y.*im_mask))];
    Covxx =  1/ellipseSize * sum(sum(((X - ellipseMean(1)).^2 ).*im_mask));
    Covyy =  1/ellipseSize * sum(sum(((Y - ellipseMean(2)).^2 ).*im_mask));
    Covxy =  1/ellipseSize * sum(sum(((X - ellipseMean(1)).*(Y - ellipseMean(2))).*im_mask));
    ellipseCov = [Covxx, Covxy; Covxy, Covyy];
    [~,D] = eig(ellipseCov);
    a=2 * sqrt(diag(D));
    ellipseMeanstack(i,:) = ellipseMean - [(size(projections,1)+1) / 2,(size(projections,2)+1) / 2];
    projectionDir(:,i) = 1/a(2) * [sqrt(a(2)^2 -4 * Covxx), sqrt(a(2)^2 -4 * Covyy), a(1)];
    if sign(projectionDir(1,i)*projectionDir(2,i)) ~= -sign(Covxy)
        projectionDir(1,i)=-projectionDir(1,i);
    end
    est_scales(i) = a(2);
end
projectionDir= projectionDir';
est_scales = est_scales/mean(est_scales);
end

