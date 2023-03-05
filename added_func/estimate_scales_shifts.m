function [ unscaled_centered_projections, iterative_LS_est_scales, iter ,...
    clstack, corrstack, iterative_LS_est_shifts, ref_q,ref_shifts,ref_scales] = estimate_scales_shifts(projections, max_shift_factor, max_rescale_factor,...
    mask_radius, n_r, n_theta, clusteringThresh, ref_q,ref_shifts,ref_scales)
% This function estimates the scales and unscales the projections
if exist('ref_q')
    isref_q=1;
else
    isref_q=0;
    ref_q=[];
end
if exist('ref_shifts')
    isref_shifts=1;
else
    isref_shifts=0;
    ref_shifts=[];
end
if exist('ref_scales')
    isref_scales=1;
else
    isref_scales=0;
    ref_scales=[];
end

Nscales=10;
Nshifts=10;
tol=1e-3;
maxiter=20;
iterative_LS_est_scales=ones([size(projections,3),1]);
iterative_LS_est_shifts=zeros([size(projections,3),2]);
unscaled_centered_projections=projections;
iter=0;
scale_err_crit=1;
shift_err_crit=1;

fprintf('========== Begin iterative shift/scale estimation loop: maxiter=%d, tolerance=%1.1e) ==========\n', maxiter, tol);
while( scale_err_crit>tol && shift_err_crit>0.99 && iter<maxiter)
    
    % Mask projections
    [np,~]=mask_fuzzy(unscaled_centered_projections, mask_radius);
    
    % Compute polar Fourier transform, using radial resolution n_r and angular
    % resolution n_theta. n_theta is the same as above.
    [npf,~]=cryo_pft(np,n_r,n_theta);
    
    [clstack,corrstack, shift_equations,~,scale_equations,~] = ...
        commonlines_gaussian_scale(npf,max_shift_factor,Nshifts, Nscales, max_rescale_factor);
    
    % Remove projections that aren't sufficiently similar
    COV = corrstack + corrstack' + eye(size(corrstack,1));
    
    PC = pcacov(COV);
    valid_eq = abs(median(COV,2)-1)<clusteringThresh;
    valid_eq =  abs(COV*PC(:,1)/max(COV*PC(:,1))-1)<clusteringThresh;
    %     viewstack(projections,5,5);
    %     L = eye(size(COV,1)) - bsxfun(@rdivide, COV, sum(COV,1));
    %
    %     D = diag(sum(COV,1));
    %     L1 = eye(size(COV,1)) - D^(-0.5) * COV * D^(-0.5);
    %     P = bsxfun(@rdivide, COV, sum(COV,1));
    %     [V,D] = eig(L);
    %     [~,sort_idx] = sort(diag(D));
    %     V2 = V(:,sort_idx(2));
    %     V3 = V(:,sort_idx(3));
    %     V4 = V(:,sort_idx(4));
    %     X1 = COV*V2;
    %     X2 = COV*V3;
    %     X3 = COV*V4;
    %
    %     figure; scatter(X1,X2,5,group_labels);
    %     cmap = hsv(max(group_labels));    %or build a custom color map
    %     colormap( cmap );
    %
    %     dataPts = [X1'; X2'; X3'];
    %     bandWidth = 0.25;
    %     plotFlag = false;
    %     [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(dataPts,bandWidth,0);
    %     [kmeansIdx, centers] = kmeans(dataPts', max(group_labels));
    %     sortedKmeansIdx = [];
    %     for g = 1:max(group_labels)
    %         gidx = (group_labels == g);
    %         gcenter = mean(dataPts(:,gidx),2);
    %         dist = sqrt(sum((centers-repmat(gcenter',[size(centers,1),1])).^2,2));
    %         [~,idx] = min(dist);
    %         sortedKmeansIdx = [sortedKmeansIdx; g*ones(sum(kmeansIdx==idx),1)];
    %     end
    %     figure; scatter3(X1, X2, X3, 5,data2cluster);
    %     cmap = jet(max(data2cluster));    %or build a custom color map
    %     colormap( cmap );
    %
    %     figure; scatter3(X1,X2,X3, 5,kmeansIdx);
    %     cmap = jet(3);    %or build a custom color map
    %     colormap( cmap );
    %
    %     C = confusionmat(group_labels,sortedKmeansIdx);
    %     figure; imagesc(C);
    
    if any(not(valid_eq))
        unscaled_centered_projections = unscaled_centered_projections(:,:,valid_eq);
        projections = projections(:,:,valid_eq);
        iterative_LS_est_scales = iterative_LS_est_scales(valid_eq);
        iterative_LS_est_shifts = iterative_LS_est_shifts(valid_eq,:);
        if isref_q
            ref_q = ref_q(:,valid_eq);
        end
                if isref_shifts
            ref_shifts = ref_shifts(valid_eq,:);
                end
                if isref_scales
            ref_scales = ref_scales(valid_eq);
        end
        fprintf('Removed %i projections \n',sum(not(valid_eq)));
        continue;
    end
    
    [est_shifts, est_scales] = est_shifts_scales( shift_equations, scale_equations);
    
    
    iterative_LS_est_scales = iterative_LS_est_scales.*est_scales;
    iterative_LS_est_shifts = iterative_LS_est_shifts + est_shifts;
    centered_projections = cryo_addshifts(projections,-full(iterative_LS_est_shifts));
    unscaled_centered_projections = scale_projections(centered_projections, 1.0./iterative_LS_est_scales);
    
    iter=iter+1;
    
    max_rescale_factor = max(full(est_scales));
    max_shift_factor = full(max(sqrt(sum(abs(est_shifts).^2,2))));
    scale_err_crit = abs(max_rescale_factor-1);
    fprintf('Iteration #%d: scale error criteria=%f  shift=%f \n',...
        iter, scale_err_crit, max_shift_factor);
end

end



