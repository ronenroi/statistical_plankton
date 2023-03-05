function [ clstack,correlation, shift_equations,shift_equations_map, scale_equations, scale_equations_map,shifts_1d]...
                        = commonlines_gaussian_scale_roi( pf,rotations,max_shift,n_shift, Nscales, max_scale,sigma)
% Detect the commonlines by searching for the maximun cross-correlation 
% between the rays on the polar Fourier Transforms. Gaussian filter applied 
% to damp the noise in the high frequency domain.

% Input: 
%
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%
%
% Output:
%
%   clstack     Common lines matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2. 
%   corrstack   The correlation of the common line between projections k1
%       and k2. Since corrstack is symmetric, it contain entries only above
%       the diagonal. corrstack(k1,k2) measures how ''common'' is the between
%       projections k1 and k2. Large value means high-similariry.
%   shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%   shift_equations_map   2D array of size n_proj by n_proj. Entry (k1,k2)
%       is the index of the equation (row number) in the array
%       "shift_equations" that corresponds to the common line between
%       projections k1 and k2. shift_map is non-zero only for k1<k2. 
%
% 
%
% This is a revision on cryo_clmatrix_v3.m by Yoel.
% Added feature: accelerate the searching process by comparing one slice
% with all the other slices simultaneously in each iteration. 
%
% Lanhui Wang, July 2, 2013


scales = exp(linspace(-log(max_scale), log(max_scale), Nscales));

if nargin<2
    max_shift=15; % Maximal shift between common-lines in pixels. The 
                  % shift  is from -max_shift to max_shift. 
end

if nargin<3
    n_shift=10; % Resolution of shift estimation in pixels.
end

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
                               
shift_equation_idx=1;  % The equation number we are currently processing.
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
for k1=1:n_theta
    for l=1:n_proj
        coefficients(:,k1,l)=coefficients(:,k1,l)/norm(coefficients(:,k1,l));
    end
    
end
C=coefficients(:,1:n_theta/2,:);
C=reshape(C,rmax,n_theta/2*n_proj); % stack all the Fourier slices
rk2=rk(1:rmax);

% Accelerate the searching process by comparing one slice
% with all the other shifted slices in each iteration.
for k1=1:n_proj
    correlation=-999*ones(n_proj,n_proj);
    idx = zeros([1, n_proj-k1]);
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(coefficients(:,:,k1), size(coefficients(:,:,k1)).*[M, 1]);
        if M<1.0
            scale_coeff = zeros(size(coefficients(:,:,k1)));
            scale_coeff(end-size(scaled_mag_coeff,1)+1:end,:)=scaled_mag_coeff;
        else
            scale_coeff = scaled_mag_coeff(end-N+1:end,:);
        end
         
         for kk=1:size(scale_coeff,2)
            scale_coeff(:,kk) = (1/norm(scale_coeff(:,kk)))*scale_coeff(:,kk);
        end
        temp_coef=zeros(N,n_theta,n_shift);

        % Generate all the shifted copies of k_th slice.
        for shiftidx=1:n_shift
            shift=-max_shift+(shiftidx-1)*shift_step;
            shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
            temp_coef(:,:,shiftidx)=bsxfun(@times,scale_coeff,shift_phases);
        end
        %%
    Ri=rotations(:,:,k1);
for k2=k1+1:n_proj
    Rj=rotations(:,:,k2);
    [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
    
    % To match cryo_clmatrix, cij is always less than PI and cji may be be
    % larger than PI.
     cl_i_vec = mod(cij-sigma : cij+sigma,n_theta) +1;
     cl_j_vec = mod(cji-sigma : cji+sigma,n_theta) +1;
     
    temp_coef_i = temp_coef(:,cl_i_vec,:);
     temp_coef_j = coefficients(:,cl_j_vec,k2);
     temp_coef_i=reshape(temp_coef_i,N,n_shift*length(cl_j_vec));
     corr = temp_coef_i' * temp_coef_j;
     corr=reshape(corr,n_shift*length(cl_j_vec)^2,1);
     [temp_correlation, temp_idx] = max(real(corr)); % pick the largest one.
     if temp_correlation>correlation(k1,k2)
            scale_1d(k1, k2) = log(M);
        correlation(k1,k2) = temp_correlation;
            [cl1, est_shift, cl2]=ind2sub([length(cl_j_vec) n_shift length(cl_j_vec)],temp_idx);
            cl1=cl_i_vec(cl1);
            cl2=cl_j_vec(cl2);
    cl2(cl1>n_theta/2)=cl2(cl1>n_theta/2)+n_theta/2;
    cl1(cl1>n_theta/2)=cl1(cl1>n_theta/2)-n_theta/2;
    shifts_1d(k1,k2)=-max_shift+(est_shift-1)*shift_step;
    clstack(k1,k2)=cl1;
    clstack(k2,k1)=cl2;
     end
end
%%
%         temp_coef=reshape(temp_coef,N,n_shift*n_theta);
% 
%         corr=temp_coef'* C(:,n_theta/2*k1+1:end); % Compute the cross correlations with all other slices.
%         corr=reshape(corr,n_shift*n_theta^2/2,n_proj-k1);
%         
%         [temp_correlation, temp_idx] = max(real(corr)); % pick the largest one.
%         ind1 = temp_correlation>correlation;
%         
%         scale_1d(k1, padarray(ind1, [0,k1], 'pre')) = log(M);
%         correlation(ind1) = temp_correlation(ind1);
%         idx(ind1) = temp_idx(ind1);
    end
    

end

%%
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