function [vol,RsCn,est_shifts,est_scales] = reconstruct_cn(projs,rots,n_symm,n_r,n_theta,max_shift,shift_step,max_scale,Nscales)

if ~exist('shift_step','var')
    shift_step = 0.5;
end

if ~exist('max_shift','var')
    max_shift = 15;%ceil(size(projs,1)*0.15); % max_shift is 15% of the image size
end

if ~exist('n_theta','var')
    n_theta = 360;
end

if ~exist('n_r','var')
    n_r = ceil(size(projs,1)*0.75);
end

[pf,~] = cryo_pft(projs,n_r,n_theta,'single');  % take Fourier transform of projections
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
[est_shifts,~,est_scales,~] = cryo_estimate_shifts_scales(pfCn,RsCn,max_shift,shift_step,max_scale,Nscales,10000,[],0);

%[est_shifts,~] = cryo_estimate_shifts(pfCn,RsCn,ceil(2*sqrt(2)*max_shift),shift_step,10000,[],0);
n = size(projs,1);
% projsCn = zeros([size(projs),n_symm]);
projsCn = [];
for i=1:n_symm
    projsCn = cat(3,projsCn,projs);
%     projsCn(:,:,i) = projs;
end
centered_projections = cryo_addshifts(projsCn,-est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./est_scales);

log_message('Reconstructing');
[ v1, ~, ~ ,~, ~, ~] = recon3d_firm( unscaled_centered_projections,RsCn,[], 1e-6, 200, zeros(n,n,n));
%[v1] = astra_reconstruction(unscaled_centered_projections,RsCn,[n,n,n],1000);
% params = struct();
% params.rot_matrices = RsCn;
% params.ctf = ones(n*ones(1, 2));
% params.ctf_idx = ones(size(projsCn, 3), 1);
% params.shifts = est_shifts.';
% params.ampl = ones(size(projsCn, 3), 1);
% 
% mean_est_opt.max_iter = 100;
% mean_est_opt.rel_tolerance = 1.0e-6;
% mean_est_opt.verbose = false;
% mean_est_opt.precision = 'single';
% 
% basis = dirac_basis(size(projsCn, 1)*ones(1, 3));
% 
% v1 = cryo_estimate_mean(single(projsCn), params, basis, mean_est_opt);

ii1=norm(imag(v1(:)))/norm(v1(:));
log_message('Relative norm of imaginary components = %e\n',ii1);
vol=real(v1);
vol(vol<0)=0;

end