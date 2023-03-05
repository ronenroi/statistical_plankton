clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
%%
scale_estimation = true;

dataFolder = '../../../data sets/Pyramimonas_longicauda'; projectionType = 'BrightField';  

projections = func_loadPlanktonData(dataFolder, projectionType);
%projections = permute(projectionsT,[2 1 3]);
%viewstack(projections,8,8);   % Display the proejctions.


%% Shift and possible scale estimation (according to flag)
n_theta=180; % Angular resolution - number of sinograms computed for each 
            % projection. This corresponds to a resolution of 5 degrees.


% Find common lines between projections. Note that the max_shift used is
% much larger than above. This is because it is a 1D shift used to seearch
% for common lines. If d is the maximal expected 2D shift in a projection
% (in any direction), then the 1D shift used to find common lines should
% equal ceil(2*sqrt(2)*d).  
max_shift = 15; 
shift_step = 1;
% n_r=33; 
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));


% estimate shifts and scales with iterative refinement
max_scale = 2;
clusteringThresh = 1;
nproj = size(projections,3);
[unscaled_centered_projections, iterative_LS_est_scales, ...
 niter, clstack, corrstack, est_shifts] = estimate_scales_shifts(projections, max_shift,...
                                                      max_scale, mask_radius, n_r, n_theta, clusteringThresh,zeros(nproj),zeros(nproj));
nproj = size(unscaled_centered_projections,3);



%viewstack(unscaled_projections,9,8);   % Display the proejctions.
%% Assign orientation using common lines, using least squares method.
% The resulting MSE should be small (of the order of 1e-4).
%[est_inv_rots_LS] = est_orientations_LS(clstack, n_theta);
[est_inv_rots_LUD] = est_orientations_LUD(clstack, n_theta); % Least Unsquared Deviations

%% 3D inversion

if scale_estimation==true
    img_size = size(unscaled_centered_projections,1);
    [ v, v_b, kernel ,err, iter, flag] = recon3d_firm( unscaled_centered_projections,...
        est_inv_rots_LUD, [], 1e-6, 200, zeros(img_size,img_size,img_size));
else
    img_size = size(projections,1); 
    [ v, v_b, kernel ,err, iter, flag] = recon3d_firm_parallel( projections,...
        est_inv_rots_LS, [], 1e-6, 200, zeros(img_size,img_size,img_size));
end

v=real(v); v(v<0)=0;
%%
save('vol_est_star_processed', 'est_inv_rots_LUD', 'unscaled_centered_projections', 'v','iterative_LS_est_scales', 'est_shifts', 'projections')
%% Refinement
if 0
    params.N = 1;
    params.defidx = ones([nproj,1]);
    params.max_shifts = 7;   % maximum shift search range
    params.c = 1;
    filename = 'refinement/refined_model';
    iter_max = 10;

    tol = 0.01;
    CTF_flag = 1;

    [ v_new ] = Refine(v, unscaled_centered_projections, params, iter_max, tol, CTF_flag, filename );
    v_new(v_new<0)=0;
end
%% Save results
save('POP_Pyramimonas_longicauda_noprune')


