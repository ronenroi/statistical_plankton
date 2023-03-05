

%% new method
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
%%
%% Load projections

load('simulation_est/simulation_estimation_LUD_scale1.mat', 'data')
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = data.max_shift;
n_theta=90;
K=size(projections,3);

mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));

%% Orientation assigment
% Use sync3N
VOTING_TICS_WIDTH=1;
J_EIGS=4;
J_WEIGHTS=true;
S_WEIGHTS=true;
Nscales=10;
Nshifts=10;
shift_step = ceil(max_shift/Nshifts);
log_message('Saving common lines');
    [masked_projs,~]=mask_fuzzy(projections,mask_radius);        %% Mask projections
    log_message('Computeing polar Fourier transform of projections. n_theta=%d, n_r=%d',n_theta,n_r);
    
    ref_clmatrix=clmatrix_cheat(q_to_rot(data.ref_q),360);

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
%%    
xbins3 = -90:1:90;
 [thetadiff ]=comparecl_model( clstack, ...
                ref_clmatrix, n_theta );
            thetadiff(logical(eye(size(thetadiff)))) = []; % Or A = A(~eye(size(A)))
            thetadiff=thetadiff-360/n_theta+1;
    N_lines = length(thetadiff);
            
            f = @(par,x)(N_lines * par(1))  * 1/(par(2)*sqrt(2*pi)) * (exp(-x.^2 ./(2*par(2)^2))+exp(-(180-x).^2./(2*par(2)^2))) ...
                + (N_lines * (1- par(1))) / 180;
            B0 = [0.8, 1.4];
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'cauchy';
            h=hist(thetadiff,xbins3);
            B = nlinfit(xbins3(1:180), h(1:180), f, B0);
            figure;
            hist(thetadiff,xbins3)
            hold on
            plot(xbins3, f(B,xbins3), '-r', 'LineWidth',1.5)
            hold off
            legend(sprintf('P %f sigma %f',B(1),B(2)));
            P=B(1);
            sigma = B(2);
            
            
            %%
            Nequations = ceil(K*(K-1)/2);
            [pairsI,pairsJ]=meshgrid(1:K,1:K);
idxI=pairsI(pairsJ>pairsI);
idxJ=pairsJ(pairsJ>pairsI);
Icl=[idxI,idxJ];
for shift_equation_idx=1:Nequations
    % Find the pair of projections participating in the current common line
    % pair.
    idxi=Icl(shift_equation_idx,1);  % Index of projection i in the pair.
    idxj=Icl(shift_equation_idx,2);  % Index of projection j in the pair.
%W_vec(shift_equation_idx) = W(idxi,idxj)/n_projs;
    % Extract the indices of the common line between Pi and Pj.
    Ri=rotations(:,:,idxi);
    Rj=rotations(:,:,idxj);
    [cij,cji]=commonline_R(Ri.',Rj.',360);
    clstack_refined(idxi,idxj)=cij;
    clstack_refined(idxj,idxi)=cji;
end
%%
            xbins3 = -90:1:90;
            [thetadiff ]=comparecl_model( clstack_refined, ...
                ref_clmatrix, 360 );
            thetadiff(logical(eye(size(thetadiff)))) = []; % Or A = A(~eye(size(A)))
             thetadiff=thetadiff+1;
            N_lines = length(thetadiff);
            
            f = @(par,x)(N_lines * par(1))  * 1/(par(2)*sqrt(2*pi)) * (exp(-x.^2 ./(2*par(2)^2))+exp(-(180-x).^2./(2*par(2)^2))) ...
                + (N_lines * (1- par(1))) / 180;
            B0 = [0.8, 1.4];
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'cauchy';
            h=hist(thetadiff,xbins3);
            B = nlinfit(xbins3(1:180), h(1:180), f, B0);
            figure;
            hist(thetadiff,xbins3)
            hold on
            plot(xbins3, f(B,xbins3), '-r', 'LineWidth',1.5)
            hold off
            legend(sprintf('P %f sigma %f',B(1),B(2)));
            P_refine=B(1);
            sigma_refine = B(2);
%%
[Rij0, r_valid_ks, r_good_ks, peak_width] = cryo_sync3n_estimate_all_Rijs...
        (clstack_refined, n_theta, VOTING_TICS_WIDTH);
    log_message('Stage 1 done: relative rotations');
    
    
    % J-synchronization
    verbose = 2;
    [J_sync,J_significance,J_eigs,J_sync_iterations,~] = ...
        cryo_sync3n_Jsync_power_method (Rij0, J_EIGS, J_WEIGHTS, verbose);
    Rij = cryo_sync3n_flip_handedness(J_sync, Rij0);
    log_message('Stage 2 done: J-synchronization');
        
        % Estimate weights for the 3x3 blocks of S
        [W_refined, Pij, scores_hist_refined] = cryo_sync3n_syncmatrix_weights(Rij0);

%%
    [pf_fine,~]=cryo_pft(masked_projs,n_r,360,'single');  % take Fourier transform of projections
    %%
[clstack_final,corrstack,shift_equations,shift_equations_map,clstack_mask] = cryo_clmatrix_cpu_roi(pf_fine,max_shift,shift_step,clstack_refined,ceil(sigma));
%%
            xbins3 = -90:1:90;
            [thetadiff ]=comparecl_model( clstack_final, ...
                ref_clmatrix, 360 );
            thetadiff(logical(eye(size(thetadiff)))) = []; % Or A = A(~eye(size(A)))
             thetadiff=thetadiff;
            N_lines = length(thetadiff);
            
            f = @(par,x)(N_lines * par(1))  * 1/(par(2)*sqrt(2*pi)) * (exp(-x.^2 ./(2*par(2)^2))+exp(-(180-x).^2./(2*par(2)^2))) ...
                + (N_lines * (1- par(1))) / 180;
            B0 = [0.8, 1.4];
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'cauchy';
            h=hist(thetadiff,xbins3);
            B = nlinfit(xbins3(1:180), h(1:180), f, B0);
            figure;
            hist(thetadiff,xbins3)
            hold on
            plot(xbins3, f(B,xbins3), '-r', 'LineWidth',1.5)
            hold off
            legend(sprintf('P %f sigma %f',B(1),B(2)));
            P_refine=B(1);
            sigma_refine = B(2);
%%
    log_message('Estimating shifts and scales');
    [est_shifts,~,est_scales,~,refined_corrstack]=SEM_estimate_shifts_scales_weighted(pf,rotations,max_shift,shift_step,max_scale,Nscales,ones(size(A)),10000,[],0);

    log_message('Finished estimating shifts and scales');
    iterative_LS_est_scales = iterative_LS_est_scales.*est_scales;
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    
    centered_projections = cryo_addshifts(projections,-iterative_LS_est_shifts);
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./iterative_LS_est_scales);
    
    
    
    max_scale = max(abs(1-est_scales))+1;
%     max_shift = max(sqrt(sum(abs(est_shifts).^2,2))) ;
    
    %     max_scale = max_scale/1.1;
    max_shift = max_shift - 1.5;
    scale_err_crit = abs(max_scale-1);
    fprintf('Iteration #%d: scale error criteria=%f  shift=%f \n',...
        iter, scale_err_crit, max_shift);
    if 0
        [MSE(iter), err, O, aligned_rots] = check_MSE(rotations,data.ref_q);
        MSE_shifts(iter) = mean(sqrt(sum(abs(iterative_LS_est_shifts - data.ref_shifts).^2,2))) / 2;
        MSE_scales(iter) = mean(abs(1-(iterative_LS_est_scales(:) ./ data.ref_scales(:))));
        max_shifts(iter) = max(sqrt(sum(abs(iterative_LS_est_shifts - data.ref_shifts).^2,2))) / 2;
        Max_scales(iter) = max(abs(1-(iterative_LS_est_scales(:) ./ data.ref_scales(:))));
    end
    save(outparams,'-append');
    




