function [est_shifts,est_scales] = max_estimate_shifts_scales_cn(pf,rots,n_symm,max_shift,shift_step,max_scale,Nscales)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 15;%ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end


pfCn = [];
for i=1:n_symm
    pfCn = cat(3,pfCn,pf);
end

g = [cosd(360/n_symm) -sind(360/n_symm) 0; 
	 sind(360/n_symm)  cosd(360/n_symm) 0; 
	 0 				 0  1]; % a rotation of 90 degrees about the z-axis

nImages = size(rots,3);

RsCn = zeros(3,3,nImages*n_symm);

for k=1:nImages
    rot = rots(:,:,k);
    for s=0:n_symm-1
        RsCn(:,:,s*nImages+k) = g^s*rot; 
    end
end
max_shift =ceil(2*sqrt(2)*max_shift);
log_message('Estimating shifts and scales');
[est_shifts_Cn,~,est_scales_Cn,~,corrstack] = estimate_shifts_scales_cn(pfCn,RsCn,max_shift,shift_step,max_scale,Nscales,10000,[],0);

COV = corrstack + corrstack' + eye(size(corrstack,1));
    
    PC = pcacov(COV);
    scores =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1);
    est_shifts = zeros(nImages,2);
    est_scales = zeros(nImages,1);
for i=1:nImages
    [~,ind] = max(scores(i:nImages:end));
    ind = ind - 1;
    est_shifts(i,:) = est_shifts_Cn(ind*nImages+i,:);
    est_scales(i) = est_scales_Cn(ind*nImages+i);
end

end