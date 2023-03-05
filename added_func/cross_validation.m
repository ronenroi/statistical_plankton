function [MSE,tested_projection_stack] = cross_validation(unscaled_centered_projections,inv_rot,iter)
ind=0;
for i=1:size(unscaled_centered_projections,3)
    ind =ind +1;
 cross_validation_projections = unscaled_centered_projections;
tested_projection = cross_validation_projections(:,:,i) ;
cross_validation_projections(:,:,i) = [];
rotations = inv_rot;
tested_rot = rotations(:,:,i)';
rotations(:,:,i) = [];

s = size(unscaled_centered_projections,1);
vol_size = [s s s];

[vol] = astra_reconstruction(cross_validation_projections,rotations,vol_size,iter);
vol(vol<0) = 0;
projection=cryo_project(vol,tested_rot,size(vol,1),'single'); % generate phantom projecitons
projection=permute(projection,[2 1]);   % col x n_proj x row
projection(projection<0)=0;
projection_stack(:,:,ind) = projection;
tested_projection_stack(:,:,ind) = tested_projection;
% figure;subplot(121);imagesc(projection)
% subplot(122);imagesc(tested_projection)
MSE(i) = mean((projection(:) - tested_projection(:)).^2);
end
save('3DPMPO_projection','projection_stack','tested_projection_stack','MSE')

end

