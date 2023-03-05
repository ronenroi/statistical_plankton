function [ MSE ,err_vec,recon_proj] = MSEvol( vol,unscaled_centered_projections,inv_rot )

%% projections comparison
n=size(unscaled_centered_projections,3);                              % projection number
rot=zeros(size(inv_rot));
MSE=zeros(1,n);
for i=1:n
    rot(:,:,i)=inv_rot(:,:,i)';
end
recon_proj = double( cryo_project(vol,rot,size(vol,1),'single')); % generate projecitons
recon_proj = permute(recon_proj,[2,1,3]);
recon_proj(recon_proj<0)=0;
for i=1:n
    m=unscaled_centered_projections(:,:,i);
    s = recon_proj(:,:,i);
  
t= mean((s(:) - m(:)).^2);


err_vec(i) = (max(max( unscaled_centered_projections(:,:,i))));


end

end

