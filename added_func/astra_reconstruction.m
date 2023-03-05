function [vol] = astra_reconstruction(projections,rotations,vol_size,iter)

% if simulation
% projsCn=permute(unscaled_centered_projections,[1 3 2]);   % col x n_proj x row
% if real data
projections=permute(projections,[2 3 1]);   % col x n_proj x row
vectors = rot2astraVec(rotations);
vol_geom = astra_create_vol_geom(permute(vol_size,[2,3,1]));
proj_geom = astra_create_proj_geom('parallel3d_vec', size(projections,3), size(projections,1), vectors);
sino_id = astra_mex_data3d('create','-proj3d', proj_geom, projections);
% store volume
vol_id = astra_mex_data3d('create','-vol', vol_geom, 0);
% create sinogram
cfg = astra_struct('SIRT3D_CUDA');
cfg.ProjectionDataId = sino_id;
cfg.ReconstructionDataId = vol_id;
cfg.option.MinConstraint = 0;
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('iterate', alg_id, iter);
astra_mex_algorithm('iterate', alg_id);
astra_mex_algorithm('delete', alg_id);
vol = astra_mex_data3d('get',vol_id);
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', vol_id);

end

