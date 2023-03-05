function outProjections = scale_projections ( projections, scales )
% This function function unscales a stack of projections according to their
% absolute scales and resizes them to fit a single image size

projections(projections<0)=0;
[imsize,~,nproj] = size(projections);
Rin = imref2d([imsize, imsize]);
Rin.XWorldLimits = Rin.XWorldLimits-mean(Rin.XWorldLimits);
Rin.YWorldLimits = Rin.YWorldLimits-mean(Rin.YWorldLimits);

outProjections = zeros([imsize,imsize,nproj]);
gcp;
parfor ii=1:nproj

    tform = affine2d([scales(ii)  0 0; 
                     0 scales(ii) 0; 
                     0 0 1]);
    outProjections(:,:,ii) = imwarp(projections(:,:,ii), Rin, tform, 'OutputView',Rin);
end

end


