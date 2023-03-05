%% reco with gt values
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../');
initpath;
cd(myPath);
%%
load('simulation_est/simulation_estimation_LUD_scale1.mat', 'data')
[MSE, err, O, aligned_rots] = check_MSE(q_to_rot(data.ref_q),data.ref_q);
centered_projections = cryo_addshifts(data.projections,-data.ref_shifts);
unscaled_centered_projections = scale_projections(centered_projections, 1.0./(data.ref_scales));
img_size = size(unscaled_centered_projections,1);
[ v1, v_b, kernel ,err, iter, flag] = recon3d_firm( unscaled_centered_projections,...
    data.ref_inv_rot_matrices,[], 1e-6, 200, zeros(img_size,img_size,img_size));  
est_vol_new=real(v1);
est_vol_new(est_vol_new<0)=0;
 v=est_vol_new;
    beta = vol;
    padlen = (size(v,1) - size(beta,1))/2;
    p = padarray(beta,[padlen,padlen,padlen]);
    corr = fftshift(ifftn(fftn(ifftshift(v)).*conj(fftn(ifftshift(p))),'symmetric'));
    [~,midx]=max(corr(:));
    [x1,y1,z1] = ind2sub(size(p),midx);
    x2 = size(beta) - 1;
    eps = inf;
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
                    est_vol_new_register = v_new_registered;
                end
                eps = min(eps, new_eps);
            end
        end
    end
    
    save('gt_recon','est_vol_new_register','eps');
    
