%% Basic abinitio reconstruction - example 1.
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
register=1;
%% Load and display projections
load cleanrib.mat
vol = volref;
max_scales = linspace(1.0,6,20);
max_scales = max_scales(5); %take only max scale of 2

%%
for kk = 1:length(max_scales)

max_scale = max_scales(kk);
data = load(['simulated_data\scale' num2str(2) '_rotations_13_11.mat']);
projections = data.projections;
beta = vol;
scale_std(kk) = std(data.ref_scales);

%  viewstack(projections,5,5);   % Display the proejctions.

%% Compute common lines from projections

n_theta=180; % Angular resolution - number of sinograms computed for each 
            % projection. This corresponds to a resolution of 2 degrees.

% Find common lines between projections. Note that the max_shift used is
% much larger than above. This is because it is a 1D shift used to seearch
% for common lines. If d is the maximal expected 2D shift in a projection
% (in any direction), then the 1D shift used to find common lines should
% equal ceil(2*sqrt(2)*d).  
max_shift = ceil(2*sqrt(2)*data.max_shift) ;
shift_step = 1;
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));

% [np,~]=mask_fuzzy(projections, mask_radius);
% [npf,~]=cryo_pft(np,n_r,n_theta);
% [ clstack_SPR,~, shift_equations,~] = ...
%         commonlines_gaussian(npf,max_shift,shift_step ); %w.o scales and shifts estimation

% estimate shifts and scales with iterative refinement
clusteringThresh = 1;


[unscaled_centered_projections, iterative_LS_est_scales, ...
 ~, clstack_GSPR, ~, iterative_LS_est_shifts, ref_q_GSPR,valid_shifts,valid_scales] = estimate_scales_shifts(projections, max_shift,...
                                                      max_scale, mask_radius, n_r, n_theta, clusteringThresh, data.ref_q,data.ref_shifts,data.ref_scales);
%[est_inv_rots_LS_GSPR] = est_orientations_LS(clstack_GSPR, n_theta);
% [est_inv_rots_LUD_SPR] = est_orientations_LUD(clstack_SPR, n_theta); % Least Unsquared Deviations
[est_inv_rots_LUD_GSPR] = est_orientations_LUD(clstack_GSPR, n_theta); % Least Unsquared Deviations

%[MSE_LS, err_LS, O1] = check_MSE(est_inv_rots_LS, data.ref_q);
% [MSE_LUD_SPR(kk), ~, ~, aligned_rots_SPR] = check_MSE(est_inv_rots_LUD_SPR, data.ref_q);
[MSE_LUD_GSPR(kk), ~, ~, aligned_rots_GSPR_LUD] = check_MSE(est_inv_rots_LUD_GSPR, ref_q_GSPR);
%[MSE_LS_GSPR(kk), ~, ~, aligned_rots_GSPR_LS] = check_MSE(est_inv_rots_LS_GSPR, ref_q_GSPR);

% disp(['MSE LUD SPR: ' num2str(MSE_LUD_SPR(kk))])
disp(['MSE LUD GSPR: ' num2str(MSE_LUD_GSPR(kk))])
% disp(['MSE LS: ' num2str(MSE_LS_GSPR(kk))])

%%
%if MSE_LUD_GSPR(kk) < MSE_LS_GSPR(kk)
    aligned_rots_GSPR = aligned_rots_GSPR_LUD;
% else
%         aligned_rots_GSPR = aligned_rots_GSPR_LS;
% end

%%
scale_factor = iterative_LS_est_scales' * valid_scales / (iterative_LS_est_scales' * iterative_LS_est_scales);
    centered_projections = cryo_addshifts(unscaled_centered_projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
     [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);
 [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
     [est_shifts,~]=cryo_estimate_shifts(pf,aligned_rots_GSPR,max_shift,shift_step,10000,[],0);
     iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
 centered_projections = cryo_addshifts(unscaled_centered_projections,-iterative_LS_est_shifts);
     unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
%    
%%
% img_size_SPR = size(projections,1);
img_size_GSPR = size(unscaled_centered_projections,1);

    [ v_GSPR, v_b, kernel ,err, iter, flag] = recon3d_firm( unscaled_centered_projections,...
aligned_rots_GSPR,[], 1e-6, 200, zeros(img_size_GSPR,img_size_GSPR,img_size_GSPR));
%     
%     
%  [ v_SPR, v_b, kernel ,err, iter, flag] = recon3d_firm_parallel( projections,...
%          aligned_rots_SPR, [], 1e-6, 200, zeros(img_size_SPR,img_size_SPR,img_size_SPR));

v_GSPR=real(v_GSPR); v_GSPR(v_GSPR<0)=0;
% v_SPR=real(v_SPR); v_SPR(v_SPR<0)=0;

%% Save results
save(['simulation_est\simulation_estimation_LUD_scale' num2str(kk)], 'unscaled_centered_projections')
save(['simulation_est\simulation_estimation_LUD_scale' num2str(kk)],...
    'data', 'v_GSPR', 'aligned_rots_GSPR', 'clusteringThresh','clstack_GSPR','ref_q_GSPR','MSE_LUD_GSPR')%, 'v_SPR','aligned_rots_SPR','clstack_SPR','MSE_LUD_SPR')
    save(['simulation_est\simulation_estimation_LUD_scale' num2str(kk)],'scale_factor','-append');




%% Register reconstruction
if register
recoveries{1} = v_GSPR; 
%recoveries{2} = v_SPR; 

for ii=1:2
    v=recoveries{ii};
    padlen = (size(v,1) - size(beta,1))/2;
    p = padarray(beta,[padlen,padlen,padlen]);
    corr = fftshift(ifftn(fftn(ifftshift(v)).*conj(fftn(ifftshift(p))),'symmetric'));
    [~,midx]=max(corr(:));
    [x1,y1,z1] = ind2sub(size(p),midx);
    x2 = size(beta) - 1;
    eps = 999;
    for i=linspace(1,1.1,10)
        v=imresize3(recoveries{ii},i);

    for s1=-10:10
        for s2=-10:10
            for s3=-10:10
                xs = round([x1,y1,z1]-x2/2)+[s1,s2,s3];
                xe = round([x1,y1,z1]+x2/2)+[s1,s2,s3];
                if min([xs,xe])<1 || max([xs,xe])>size(v,1)
                    continue;
                end
                v_new_registered = v(xs(1):xe(1),xs(2):xe(2),xs(3):xe(3));

                new_eps = (1/sum(beta(:))) * sum(abs(beta(:)-v_new_registered(:)));
                if new_eps<eps
                    Xs{ii} = xs; Xe{ii}=xe;
                    recovery_reg{ii} = v_new_registered;
                    scale_i = i;
                end   
                eps = min(eps, new_eps);
            end
        end
    end
    
    if ii==1
        curr_eps_GSPR = eps;
        eps_GSPR(kk) = eps;
        v_GSPR_registered=recovery_reg{ii};
        fprintf('Iter: %d  eps_GSPR=%2.3f  std=%2.3f \n', kk, eps, scale_std(kk));
    else
        curr_eps_SPR = eps;
        eps_SPR(kk) = eps;
        v_SPR_registered=recovery_reg{ii};
        
        fprintf('Iter: %d  eps_SPR=%2.3f  std=%2.3f \n', kk, eps, scale_std(kk));

    end
    end

end
            save(['simulation_est\simulation_estimation_registered_LUD_scale' num2str(kk)],...
    'data', 'v_GSPR_registered', 'curr_eps_GSPR')%,'v_SPR_registered', 'curr_eps_SPR')
end
end
%%
% 
% if scale_estimation
%     eps_GSPR(kk) = eps;
% else
%     eps_SPR(kk) = eps;
% end
% 
% 
% 
% 
% %% plot error
% 
% figure; plot(scale_std, eps_GSPR, 'b'); hold on; 
% plot(scale_std, eps_SPR,'--r');
% legend('GSPR', 'SPR', 'Location', 'northwest');
% 
% %% Refinement
% 
% params.N = 1;
% params.defidx = ones([nproj,1]);
% params.max_shifts = ceil(max(est_shifts(:)));   % maximum shift search range
% params.c = 1;
% filename = 'refinement/refined_model';
% iter_max = 20;
% tol = 0.01;
% CTF_flag = 1;
% 
% [ v_new ] = Refine(v, unscaled_centered_projections, params, iter_max, tol, CTF_flag, filename );
% v_new=real(v_new); v_new(v_new<0)=0;
% 
% %% Register recovered volume to ground truth
% padlen = size(v_new,1) - size(data.ref_vols{1},1);
% p = padarray(data.ref_vols{1},[padlen/2,padlen/2,padlen/2]);
% %p = padarray(p,[1,1,1],'pre');
% corr = fftshift(ifftn(fftn(ifftshift(v_new)).*conj(fftn(ifftshift(p))),'symmetric'));
% [max_corr,midx]=max(corr(:));
% [x1,y1,z1] = ind2sub(size(p),midx);
% x2 = size(data.ref_vols{1}) - 1;
% xs = round([x1,y1,z1]-x2/2);
% xe = round([x1,y1,z1]+x2/2);
% v_new_registered = v_new(xs(1):xe(1),xs(2):xe(2),xs(3):xe(3));
% 
% %figure; vol3d('cdata',v_new_registered); figure; vol3d('cdata',data.ref_vol);
% %% Save Results
% gt = data.ref_vol;
% save('Rebuttal_cloud',...
%     'gt', 'aligned_rots_GSPR' , 'aligned_rots_SPR','unscaled_projections', 'v_GSPR', 'v_SPR', ...
%     'iterative_LS_est_scales', 'xs_GSPR','xe_GSPR', 'xs_SPR', 'xe_SPR', 'v_GSPR_Reg', 'v_SPR_Reg')
% 
% %%
% save(['Hydra_Frag488_withnoise_new' num2str(length(data.ref_vols))],...
%     'data', 'est_inv_rots_LS', 'est_inv_rots_LUD' , 'unscaled_projections', 'v', ...
%     'v_new', 'MSE_LUD','iterative_LS_est_scales', 'est_shifts','O', 'v_new_registered', 'xs','xe')
% 
% %WriteMRC(v,1,'example1.mrc'); % Output density map reconstructed from projections.

