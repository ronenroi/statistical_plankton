function SEM_worker(projections,outparams,showfigs,...
    n_theta,n_r,max_shift,shift_step,max_scale,data)
% CRYO_ABINITO_C1_WORKER  Worker function for C1 abainitio reconstruction
%
% Internal function called by abinitio reconstruction functions.
% Do not call this function directly.
%
% Executes the C1 abinitio reconstruction algorithm specified by 'algo'.
% Takes as input an MRC stack of images and reconstructs an abinito model
% from the images. The reconstructed volume is stored in outvol.
% For a detailed description of the parameters see below.
%
% Parameters

%   instack     Name of MRC file containing the projections (or class
%               averages) from which to estimate an abinitio model.
%   outvol      Name of MRC file into which to save the reconstructed
%               volume.
%   outparams   Name of MAT file in which intermediate outcomes of the
%               reconstruction algorithm are save (such as estimated
%               rotations and shifts). Used to debugging and detailed
%               analysis of the results.
%   showfigs    (optional) Set to 1 to display quaility assessment figures.
%               Currently only consistency histogram for the estimated
%               rotations is display. Default is 0.
%   verbose     (Optional) Set to 1 to print verbose long messages. Default
%               is 1.
%   ntheta      (Optional) Angular resolution for common lines detection.
%               Default 360.
%   n_r         (Optional) Radial resolution for common line detection.
%               Default is half the width of the images.
%   max_shift   (Optional) Maximal 1d shift (in pixels) to search between
%               common-lines. Default is 15% of image width of the images.
%   shift_step  (Optional) Resolution of shift estimation in pixels. Note
%               that shift_step can be any positive real number. Default:1.
%



%% Load projections
K=size(projections,3);
log_message('projections loaded. Using K=%d projections of size %d x %d',K,size(projections,1),size(projections,2));
if size(projections,1)~=size(projections,2)
    error('Input images must be square');
end



%% Orientation assigment
% Use sync3N
VOTING_TICS_WIDTH=1;
J_EIGS=4;
J_WEIGHTS=true;
S_WEIGHTS=true;
Nscales=10;
Nshifts=10;
tol=1e-3;
maxiter=20;
iterative_LS_est_scales=ones([size(projections,3),1]);
iterative_LS_est_shifts=zeros([size(projections,3),2]);
unscaled_centered_projections=projections;
iter=0;
scale_err_crit=1;
errhist = zeros(size(projections,3),maxiter);
mask_radius=round(size(projections,1)*0.75);
discard_projections=[];
log_message('Masking projections. Masking radius is %d pixels',mask_radius);
fprintf('========== Begin iterative shift/scale estimation loop: maxiter=%d, tolerance=%1.1e) ==========\n', maxiter, tol);

log_message('Saving common lines');
save(outparams,'n_theta','n_r','max_shift','shift_step','max_scale');
while( scale_err_crit>tol && max_shift>0.99 && iter<maxiter)
    iter=iter+1;
    [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);        %% Mask projections
    log_message('Computeing polar Fourier transform of projections. n_theta=%d, n_r=%d',n_theta,n_r);
    
    
    [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
    [clstack,corrstack] = only_commonlines_gaussian_shift_scale(pf,max_shift,Nshifts, Nscales, max_scale);%% Find only common lines from projections
    
    %% Estimate relative rotations
    [Rij0, r_valid_ks, r_good_ks, peak_width] = cryo_sync3n_estimate_all_Rijs...
        (clstack, n_theta, VOTING_TICS_WIDTH);
    log_message('Stage 1 done: relative rotations');
    
    
    % J-synchronization
    verbose = 2;
    [J_sync,J_significance,J_eigs,J_sync_iterations,~] = ...
        cryo_sync3n_Jsync_power_method (Rij0, J_EIGS, J_WEIGHTS, verbose);
    Rij = cryo_sync3n_flip_handedness(J_sync, Rij0);
    log_message('Stage 2 done: J-synchronization');
    
    
    % Build 3NX3N Matrix S
    S = cryo_sync3n_syncmatrix(Rij);
    log_message('Stage 3 done: constructing matrix S');
    
    
    % S Weighting
    if S_WEIGHTS
        % Estimate weights for the 3x3 blocks of S
        [W, Pij, scores_hist] = cryo_sync3n_syncmatrix_weights(Rij0);
    else
        W = ones(size(projections,3)); % Default weights are all one (no weights)
        Pij = [];
        scores_hist = struct();
    end
    log_message('Stage 4 done: computing weights');
    
    % Estimate rotations from S
    [rotations, S_eigs, ~] = cryo_sync3n_S_to_rot (S,10,W);
    
    log_message('Stage 5 done: estimating rotations');
    
    d_str=sprintf('%7.2f ',S_eigs);
    log_message('Top 10 eigenvalues of (weighted) sync matrix are %s',d_str);
    
    
    if showfigs
        clerr=cryo_syncconsistency(rotations,clstack,n_theta);
        hist(clerr(:,3),360);
        errhist(1:length(clerr(:,3)),iter) = clerr(:,3);
    end
    indeces = clerr(:,1:2);
    errors = commonline_errors(clerr(:,3),indeces,size(projections,3));
%     thr = 180;
    clusteringThresh = 1;
    COV = corrstack + corrstack' + eye(size(corrstack,1));
    
    PC = pcacov(COV);
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1)<clusteringThresh;
    %%
    A=sparse(indeces(:,1),indeces(:,2),clerr(:,3),size(pf,3),size(pf,3));
    A=full(A);
    A=A(valid_eq,:);
    A=A(:,valid_eq);
    A=A';
    A=A(logical(tril(ones(size(A')),-1))) ;
    discard_projections = cat(3,discard_projections,projections(:,:,not(valid_eq)));
    pf = pf(:,:,valid_eq);
    rotations = rotations(:,:,valid_eq);
    iterative_LS_est_scales = iterative_LS_est_scales(valid_eq);
    iterative_LS_est_shifts = iterative_LS_est_shifts(valid_eq,:);
    projections = projections(:,:,valid_eq);
    
%    
% [B,I] = sort(errors);
% figure;viewstack(unscaled_centered_projections(:,:,I),4,4)
% stack_img.(['p' num2str(iter)]) = unscaled_centered_projections;
    %%
    log_message('Num of projections after pruning %d',sum(size(pf,3)));

    
    log_message('Estimating shifts and scales');
    [est_shifts,~,est_scales,~,refined_corrstack]=SEM_estimate_shifts_scales_weighted(pf,rotations,max_shift,shift_step,max_scale,Nscales,ones(size(A)),10000,[],0);

    log_message('Finished estimating shifts and scales');
    iterative_LS_est_scales = iterative_LS_est_scales.*est_scales;
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./iterative_LS_est_scales);
    
    
    
    max_scale = max(abs(1-est_scales))+1;
    max_shift = max(sqrt(sum(abs(est_shifts).^2,2))) /2;
    
    %     max_scale = max_scale/1.1;
    %max_shift = max_shift /2;
    scale_err_crit = abs(max_scale-1);
    fprintf('Iteration #%d: scale error criteria=%f  shift=%f \n',...
        iter, scale_err_crit, max_shift);
    if 0
        [MSE(iter), err, O, aligned_rots] = check_MSE(rotations,data.ref_q);
        MSE_shifts(iter) = mean(sqrt(sum(abs(iterative_LS_est_shifts - data.ref_shifts).^2,2)));
        MSE_scales(iter) = mean(abs(1-(iterative_LS_est_scales(:) ./ data.ref_scales(:))));
        max_shifts(iter) = max(sqrt(sum(abs(iterative_LS_est_shifts - data.ref_shifts).^2,2))) ;
        Max_scales(iter) = max(abs(1-(iterative_LS_est_scales(:) ./ data.ref_scales(:))));
    end
    save(outparams,'-append');
    
end
%% global scale
if 1
    scale_factor = iterative_LS_est_scales' * data.ref_scales / (iterative_LS_est_scales' * iterative_LS_est_scales);
    save(outparams,'scale_factor','-append');
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
    [masked_projs,~]=mask_fuzzy(unscaled_centered_projections,mask_radius);
    [pf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');  % take Fourier transform of projections
    [est_shifts,~]=cryo_estimate_shifts(pf,rotations,15,shift_step,10000,[],0);
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./(iterative_LS_est_scales*scale_factor));
else
    aligned_rots = rotations ;
end

save(outparams,'-append');

%% Reconstruct downsampled volume with no CTF correction

% img_size = size(unscaled_centered_projections,1);
% [ v1, v_b, kernel ,err, iter, flag] = recon3d_firm( unscaled_centered_projections,...
%     aligned_rots,[], 1e-6, 200, zeros(img_size,img_size,img_size));
% ii1=norm(imag(v1(:)))/norm(v1(:));
% log_message('Relative norm of imaginary components = %e\n',ii1);
% est_vol_new=real(v1);
% est_vol_new(est_vol_new<0)=0;
% save(outparams,'est_vol_new','-append');
% if 0
%     [estR,estdx,estscale,est_vol_register,reflect ]=registeration(data.ref_vol,est_vol_new);
%     save(outparams,'est_vol_register','estscale','-append');
% 
% end

end
