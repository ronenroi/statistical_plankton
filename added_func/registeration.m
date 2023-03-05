function [estR,estdx,est_vol_register,reflect,estscale,eps ] = registeration( gt_vol,est_vol )
%% Register reconstruction
estR=0;
estdx=0;
reflect=0;
    eps = inf;
    est_vol_original=est_vol;
for s=1%linspace(0.,1.3,10)
    est_vol = imresize3(est_vol_original,s);
        padlen = round((size(est_vol,1) - size(gt_vol,1))/2);
    p = padarray(gt_vol,[padlen,padlen,padlen]);
    if size(p) ~= size(est_vol)
        p=p(1:size(est_vol,1),1:size(est_vol,2),1:size(est_vol,3));
    end
    corr = fftshift(ifftn(fftn(ifftshift(est_vol)).*conj(fftn(ifftshift(p))),'symmetric'));
    [~,midx]=max(corr(:));
    [x1,y1,z1] = ind2sub(size(p),midx);
    x2 = size(gt_vol) - 1;
    for s1=-10:10
        for s2=-10:10
            for s3=-10:10
                xs = round([x1,y1,z1]-x2/2)+[s1,s2,s3];
                xe = round([x1,y1,z1]+x2/2)+[s1,s2,s3];
                if min([xs,xe])<1 || max([xs,xe])>size(est_vol,1)
                    continue;
                end
                v_new_registered = est_vol(xs(1):xe(1),xs(2):xe(2),xs(3):xe(3));
                 new_eps = norm(v_new_registered(:)-gt_vol(:))/norm(gt_vol(:));
               % new_eps = (1/sum(abs(gt_vol(:)))) * sum(abs(gt_vol(:)-v_new_registered(:)));
                if new_eps<eps
                    est_vol_register = v_new_registered;
                    estscale = s;
                end
                eps = min(eps, new_eps);
            end
        end
    end
end

 %   [estR,estdx,est_vol_register,reflect]=cryo_align_densities(gt_vol,est_vol_register);

end

