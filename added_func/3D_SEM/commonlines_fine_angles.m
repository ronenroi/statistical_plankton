function [ clstack,est_shifts,est_scales]...
                        = commonlines_fine_angles( pf,rotations,max_shift,n_shifts, Nscales, max_scale,sigma)

sigma=ceil(sigma);
n_theta=size(pf,2);
n_theta2=n_theta/2;
pf=[flipdim(pf(2:end,n_theta/2+1:end,:),1) ; pf(:,1:n_theta/2,:) ];
n_projs=size(rotations,3);

Nequations = ceil(n_projs*(n_projs-1)/2); % Number of equations that will be used to estimation the shifts


 


if Nequations<n_projs
    warning('Too few equations. Increase memoryfactor. Setting Nequations to n_projs');
    Nequations=n_projs;
end

if Nequations<2*n_projs
    warning('Number of equations is small. Consider increase memoryfactor.');
end


% Allocate storage for the equations of determining the 2D shift of each
% projection. The shift equations are represented using a sparse matrix,
% since each row in the system contains four non-zeros (as it involves
% exactly four unknowns). 
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index 
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*Nequations,1);  % Row index for sparse equations system.
shift_J=zeros(4*Nequations,1);  % Column index for sparse equations system.
shift_eq=zeros(4*Nequations,1); % The coefficients of the center estimation
                                % system ordered as a single vector.                                    
shift_b=zeros(Nequations,1);    % Right hand side of the system.


%% Allocate variables used for magnification estimation
scale_1d=zeros(n_projs,n_projs);     % Estimated 1D shift between common-lines. 
scale_I=zeros(2*Nequations,1);  % Row index for sparse equations system.

scale_J=zeros(2*Nequations,1);  % Column index for sparse equations system.

scale_eq=zeros(2*Nequations,1); % The coefficients of the center estimation
    % system ordered as a single vector.
                               
scale_b=zeros(Nequations,1);   % Right hand side of the system.

% Prepare the shift phases
shift_step=2*max_shift/(n_shifts-1); % Number of shifts to try.
rmax=(size(pf,1)-1)/2;
rk=-rmax:rmax; rk=rk(:);
rk2=rk(1:rmax);
% shift_phases=zeros(rmax,n_shifts);
% for shiftidx=1:n_shifts
%     shift=-max_shift+(shiftidx-1)*shift_step;
%     shift_phases(:,shiftidx)=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
% end
scales = exp(linspace(-log(max_scale), log(max_scale), Nscales));


H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2)); 
H=H(:);  % Filter for common-line detection.

dtheta=pi/n_theta2; % Not 2*pi/n_theta, since we divided n_theta by 2 
                    % (resulting in n_theta2) to take rays of length
                    % 2*n_r-1. 

corrstack=zeros(n_projs,n_projs);    % Correlation coefficient for each common-line.

% Generate indices of Nequations pairs of common lines.
% Generate all pairs (i,j) such that j>i.
[pairsI,pairsJ]=meshgrid(1:n_projs,1:n_projs);
idxI=pairsI(pairsJ>pairsI);
idxJ=pairsJ(pairsJ>pairsI);
Icl=[idxI,idxJ];
% Pick Nequations indices from I at random.
rp=randperm(size(Icl,1));
Icl=Icl(rp(1:Nequations),:);
%Icl=Icl((1:Nequations),:);
% Iterate over the common lines pairs in Icl and for each pair find the 1D
% relative shift between the two Fourier lines in the pair.
for shift_equation_idx=1:Nequations
    % Find the pair of projections participating in the current common line
    % pair.
    idxi=Icl(shift_equation_idx,1);  % Index of projection i in the pair.
    idxj=Icl(shift_equation_idx,2);  % Index of projection j in the pair.
%W_vec(shift_equation_idx) = W(idxi,idxj)/n_projs;
    % Extract the indices of the common line between Pi and Pj.
    Ri=rotations(:,:,idxi);
    Rj=rotations(:,:,idxj);
    [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
    
    % To match cryo_clmatrix, cij is always less than PI and cji may be be
    % larger than PI.
    if cij>=n_theta/2
        cij=cij-n_theta/2;
        cji=cji-n_theta/2;
    end
    if cji<0
        cji=cji+n_theta;
    end
    
    cij=cij+1; cji=cji+1; % Since commonline_R returns zero-based indices.
    
    % Extract the Fourier rays that correspond to the common line.
    pfi=pf(:,mod(cij-sigma:cij+sigma,n_theta/2),idxi);
    
    isPfjFlipped=zeros(2*sigma+1,1);  % Is the common line in image j in the positive
                     % direction of the ray (isPfjflipped=0) or in the
                     % negative direction (isPfjflipped=1). 
                     cji_vec =  cji-sigma:cji+sigma ;
    for i=1:length(cji_vec)
        cji = cji_vec(i);
        if cji<=0
            cji = cji+n_theta;
        end
    if cji<=n_theta2
        pfj(:,i)=pf(:,cji,idxj);
    else
        pfj(:,i)=pf(:,cji-n_theta2,idxj);
        isPfjFlipped(i)=1;
    end
    end
        % Find 1D shift between pfi and pfj.
    pfi=pfi.*H;
    pfi(rmax:rmax+2,:)=0;
    pfi=cryo_raynormalize(pfi);
    pfi=pfi(1:rmax,:);  % Take half ray plus the DC
     
    pfj=pfj.*H;
    pfj(rmax:rmax+2,:)=0;
    pfj=cryo_raynormalize(pfj);
    pfj=pfj(1:rmax,:);

    

    %% rescale the smaller index ray
    if idxi < idxj
        pf_scaled = pfi;
    else
        pf_scaled = pfj;
    end
    sval = -inf;
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(pf_scaled,  size(pf_scaled) .* [M,1]);
        if M<1.0
            scale_coeff = zeros(size(pf_scaled));
            scale_coeff(end-size(scaled_mag_coeff,1)+1:end,:)=scaled_mag_coeff;
        else
            scale_coeff = scaled_mag_coeff(end-rmax+1:end,:);
        end
         for kk=1:size(scale_coeff,2)
            scale_coeff(:,kk) = (1/norm(scale_coeff(:,kk)))*scale_coeff(:,kk);
        end         
         if idxi < idxj
            pfi_scaled = scale_coeff;
            pfj_scaled = pfj;
        else
            pfj_scaled = scale_coeff;
            pfi_scaled = pfi;
         end
%%
        pfi_flipped=conj(pfi_scaled);

        for shiftidx=1:n_shifts
            shift=-max_shift+(shiftidx-1)*shift_step;
            shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
            temp_coef(:,:,shiftidx)=bsxfun(@times,pfi_scaled,shift_phases);
            temp_coef_flip(:,:,shiftidx)=bsxfun(@times,pfi_flipped,shift_phases);
        end
%%
        
         
         
%         pfi_flipped=conj(pfi_scaled);
%         pfi_stack=bsxfun(@times,pfi_scaled,shift_phases);
%         pfi_flipped_stack=bsxfun(@times,pfi_flipped,shift_phases);
        temp_coef=reshape(temp_coef,rmax,n_shifts*size(temp_coef,2));
        
        C1=2*real(temp_coef'*pfj_scaled);
        C2=2*real(pfi_flipped_stack'*pfj_scaled);
        
        [sval1,sidx1]=max(C1(:));
        [sval2,sidx2]=max(C2(:));
        if max(sval1,sval2) > sval
            scale_1d = log(M);
            if sval1>sval2 % Rays match in the same orientation.
                sval = sval1;
                dx=-max_shift+(sidx1-1)*shift_step;
            else %Rays match in opposite orientation.
                sval = sval2;
                dx=-max_shift+(sidx2-1)*shift_step;
            end
        end
        
        
        
    end
    corrstack(idxi,idxj)=sval;
    % Create a shift equation for the projections pair (k1,k2).
    idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
    shift_alpha=(cij-1)*dtheta;  % Angle of common ray in projection i.
    shift_beta= (cji-1)*dtheta;  % Angle of common ray in projection j.
    shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
    shift_J(idx)=[2*idxi-1 2*idxi 2*idxj-1 2*idxj]; % Columns of the shift variables that correspond to the current pair (k1,k2).
    shift_b(shift_equation_idx)=dx; % Right hand side of the current equation.
    
    % Compute the coefficients of the current equation.
    if ~isPfjFlipped
        shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
    else
        shift_beta=shift_beta-pi; % In the derivation we assume that all 
                        % angles are less than PI where angles larger than
                        % PI are assigned negative orientation.
        shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
    end
    
    % Scale equations
    sidx=2*(shift_equation_idx-1)+1:2*shift_equation_idx;
    scale_I(sidx)=shift_equation_idx; % Row index to construct the sparse equations.
    scale_J(sidx)=[idxi idxj]; % Columns of the shift variables that correspond to the current pair (k1,k2).
    scale_b(shift_equation_idx)=scale_1d; % Right hand side of the current equation.
    scale_eq(sidx)=[1 -1];
end

shift_equations=sparse(shift_I,shift_J,shift_eq,Nequations,2*n_projs);
shift_equations=[shift_equations shift_b(1:shift_equation_idx)];
% est_shifts=shift_equations(:,1:end-1)\shift_equations(:,end);
% est_shifts1=lsqr(shift_equations(:,1:end-1),shift_equations(:,end),1.0e-8,100);
% est_shifts = lscov([shift_equations(:,1:end-1);eye(2*n_projs)],[shift_equations(:,end);zeros(2*n_projs,1)],[W ;ones(2*n_projs,1)] );
% est_shifts=full(transpose(reshape(est_shifts,2,n_projs)));
% est_shifts1=full(transpose(reshape(est_shifts1,2,n_projs)));

scale_equations=sparse(scale_I,scale_J,scale_eq,Nequations, n_projs); 
scale_equations=[scale_equations scale_b(1:shift_equation_idx)];
% est_scales=lsqr(scale_equations(:,1:end-1),scale_equations(:,end),1.0e-8,100);
% est_scales=full(transpose(reshape(exp(est_scales),1,n_projs)));
% est_scales=lscov(scale_equations(:,1:end-1),scale_equations(:,end),W);
% est_scales=full(transpose(reshape(exp(est_scales),1,n_projs)));
[ est_shifts, est_scales] = est_shifts_scales( shift_equations, scale_equations);

if ~isscalar(shifts_2d_ref) 
    if nnz(shift_equations(:,1:end-1))<10^7
            [~,~,V]=svd(full(shift_equations(:,1:end-1)));
            s1=reshape(shifts_2d_ref.',2*n_projs,1);
            s2=reshape(est_shifts.',2*n_projs,1);
                        s3=reshape(est_shifts1.',2*n_projs,1);

            V=V(:,1:end-3); % Null space of shift_equations.
            % Compute the difference between the true and estimated shifts in
            % the subspace that is orthogonal to the null space of
            % shift_equations.

            if norm(V.'*s1)>1.0e-12
                log_message('cryo_estimate_shifts error: %8.5e',...
                    (norm(V.'*(s1-s2))/norm(V.'*s1)));
            else
                % Print absolute error
                log_message('norm(V.''*s1) = %7.5e',norm(V.'*s2));
            end
            if norm(V.'*s1)>1.0e-12
                log_message('cryo_estimate_shifts error: %8.5e',...
                    (norm(V.'*(s1-s3))/norm(V.'*s1)));
            else
                % Print absolute error
                log_message('norm(V.''*s1) = %7.5e',norm(V.'*s3));
            end
    else
        log_message('Not comparing to reference shifts - too many equations\n');
    end
end

if verbose
    if n_projs<=100
        s=svd(full(shift_equations(:,1:end-1)));
        log_message('Singular values of the shift system of equations:');
        log_message('%d  ',fliplr(s.'));
        figure;
        bar(s(end:-1:end-19))
    end
end


