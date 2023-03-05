function [clstack, corrstack, est_shifts, est_scales,shifts_1d]=SEM_estimate_shifts_scales_like_aviad_test(pf,clstack_in,...
    max_shift,n_shift,max_scale,Nscales)% Detect the commonlines by searching for the maximun cross-correlation 



scales = exp(linspace(-log(max_scale), log(max_scale), Nscales));



shift_step = 2*max_shift/(n_shift-1);

[n_r,n_theta,n_proj]=size(pf);

% pf is of size n_rxn_theta. Convert pf into an array of size
% (2xn_r-1)xn_theta, that is, take then entire ray through the origin.
pf=[flipdim(pf(2:end,n_theta/2+1:end,:),1) ; pf(:,1:n_theta/2,:) ];
temp = pf;
pf = zeros(2*n_r - 1, n_theta, n_proj);
pf(:,1: n_theta/2, :) = temp;
pf(:, n_theta/2 + 1: end, :) = flipdim(temp, 1);
pf = reshape(pf, 2*n_r - 1, n_theta * n_proj);

 
%% Apply Gaussian filter 
rmax = n_r - 1;
rk=-rmax:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2)); 
pf=bsxfun(@times,pf,H);
pf=pf(1:rmax,:);
N = rmax;
coefficients=reshape(pf,N,n_theta,n_proj);

%% Allocate variables used for shift estimation

shifts_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines. 

% Based on the estimated common-lines, construct the equations for
% determining the 2D shift of each projection. The shift equations are
% represented using a sparse matrix, since each row in the system contains
% four non-zeros (as it involves exactly four unknowns).
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index 
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*n_proj^2,1);  % Row index for sparse equations system.

shift_J=zeros(4*n_proj^2,1);  % Column index for sparse equations system.

shift_eq=zeros(4*n_proj^2,1); % The coefficients of the center estimation
    % system ordered as a single vector.
     
shift_equations_map=zeros(n_proj); % Entry (k1,k2) is the index of the 
    % euqation for the common line of projections k1 and k2. 
                               
shift_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system.
dtheta=2*pi/n_theta; 

%% Allocate variables used for magnification estimation

scale_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines. 
scale_I=zeros(2*n_proj^2,1);  % Row index for sparse equations system.

scale_J=zeros(2*n_proj^2,1);  % Column index for sparse equations system.

scale_eq=zeros(2*n_proj^2,1); % The coefficients of the center estimation
    % system ordered as a single vector.
                               
scale_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system.

%% Search for commonlines
clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.

% normalization
for k=1:n_theta
    for l=1:n_proj
        coefficients(:,k,l)=coefficients(:,k,l)/norm(coefficients(:,k,l));
    end
    
end
rk2=rk(1:rmax);
[pairsI,pairsJ]=meshgrid(1:n_proj,1:n_proj);
idxI=pairsI(pairsJ>pairsI);
idxJ=pairsJ(pairsJ>pairsI);
Icl=[idxI,idxJ];
% Accelerate the searching process by comparing one slice
% with all the other shifted slices in each iteration.
Nequations = ceil(n_proj*(n_proj-1)/2); % Number of equations that will be used to estimation the shifts
for shift_idx=1:Nequations
    
    idxi=Icl(shift_idx,1);  % Index of projection i in the pair.
    idxj=Icl(shift_idx,2);  % Index of projection j in the pair.
    % Extract the indices of the common line between Pi and Pj.
%     Ri=rotations(:,:,idxi);
%     Rj=rotations(:,:,idxj);
%     [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
cij = clstack_in(idxi,idxj);
cji = clstack_in(idxj,idxi);
    % To match cryo_clmatrix, cij is always less than PI and cji may be be
    % larger than PI.
%     if cij>=n_theta/2
%         cij=cij-n_theta/2;
%         cji=cji-n_theta/2;
%     end
%     if cji<0
%         cji=cji+n_theta;
%     end
    
%     cij=cij+1; cji=cji+1; % Since commonline_R returns zero-based indices.
    
    
    
    
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(coefficients(:,cij,idxi), size(coefficients(:,cij,idxi)).*[M, 1]);
        if M<1.0
            scale_coeff = zeros(size(coefficients(:,cij,idxi)));
            scale_coeff(end-size(scaled_mag_coeff,1)+1:end)=scaled_mag_coeff;
        else
            scale_coeff = scaled_mag_coeff(end-N+1:end);
        end
         
            scale_coeff = (1/norm(scale_coeff))*scale_coeff;
        
        temp_coef=zeros(N,n_shift);

        % Generate all the shifted copies of k_th slice.
        for shiftidx=1:n_shift
            shift=-max_shift+(shiftidx-1)*shift_step;
            shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
            temp_coef(:,shiftidx)=bsxfun(@times,scale_coeff,shift_phases);
        end
        
        temp_coef=reshape(temp_coef,N,n_shift);

        corr=temp_coef'* coefficients(:,cji,idxj); % Compute the cross correlations with all other slices.
        
        [temp_correlation, temp_shiftidx]= max(real(corr)); % pick the largest one.
        if temp_correlation>corrstack(idxi,idxj)
            scale_1d(idxi,idxj) = log(M);
           corrstack(idxi,idxj) = temp_correlation;
           shifts_1d(idxi,idxj)=-max_shift+(temp_shiftidx-1)*shift_step;

        end
        
    end
    
    clstack(idxi,idxj)=cij;
    clstack(idxj,idxi)=cji';
 end

%%
shift_equation_idx=1;

for k1=1:n_proj
    for k2=k1+1:n_proj
        idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
       
        cl1=clstack(k1,k2);
        cl2=clstack(k2,k1);
        shift_alpha=(cl1-1)*dtheta;  % Angle of common ray in projection 1.
        shift_beta= (cl2-1)*dtheta;  % Angle of common ray in projection 2.
        shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
        shift_J(idx)=[2*k1-1 2*k1 2*k2-1 2*k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        shift_b(shift_equation_idx)=shifts_1d(k1,k2); % Right hand side of the current equation
        
        % Compute the coefficients of the current equation.
        if shift_beta<pi-1e-13
            shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        else
            shift_beta=shift_beta-pi; % In the derivation we assume that all angles are less
            % than PI where angles larger than PI are assigned
            % nigative orientation.
            shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        end
        
        shift_equations_map(k1,k2)=shift_equation_idx;  % For each pair (k1,k2), store the index of its equation.
        
        
        % Scale equations
        sidx=2*(shift_equation_idx-1)+1:2*shift_equation_idx;
        scale_I(sidx)=shift_equation_idx; % Row index to construct the sparse equations.
        scale_J(sidx)=[k1 k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        scale_b(shift_equation_idx)=scale_1d(k1,k2); % Right hand side of the current equation.
        scale_eq(sidx)=[1 -1];
        
        shift_equation_idx=shift_equation_idx+1;
    end
end

shift_equation_idx=shift_equation_idx-1;
shift_equations=sparse(shift_I(1:4*shift_equation_idx),...
    shift_J(1:4*shift_equation_idx),shift_eq(1:4*shift_equation_idx),...
    shift_equation_idx,2*n_proj);
shift_equations=[shift_equations shift_b(1:shift_equation_idx)];

scale_equations_map=shift_equations_map;
scale_equations=sparse(scale_I(1:2*shift_equation_idx),...
    scale_J(1:2*shift_equation_idx),scale_eq(1:2*shift_equation_idx),...
    shift_equation_idx, n_proj); 
scale_equations=[scale_equations scale_b(1:shift_equation_idx)];
    [est_shifts, est_scales] = est_shifts_scales( shift_equations, scale_equations);
est_shifts = full(est_shifts);
est_scales = full(est_scales);
end