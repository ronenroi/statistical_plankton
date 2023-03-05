function [ clstack,corrstack]= only_commonlines_gaussian_shift_scale( pf,max_shift,n_shift, Nscales, max_scale)
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
     
                              
%% Search for commonlines
clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.

% normalization
for k=1:n_theta
    for l=1:n_proj
        coefficients(:,k,l)=coefficients(:,k,l)/norm(coefficients(:,k,l));
    end
    
end
C=coefficients(:,1:n_theta/2,:);
C=reshape(C,rmax,n_theta/2*n_proj); % stack all the Fourier slices
rk2=rk(1:rmax);

% Accelerate the searching process by comparing one slice
% with all the other shifted slices in each iteration.
for k=1:n_proj
    correlation=-inf;
    idx = zeros([1, n_proj-k]);
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(coefficients(:,:,k), size(coefficients(:,:,k)).*[M, 1]);
        if M<1.0
            scale_coeff = zeros(size(coefficients(:,:,k)));
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
        
        temp_coef=reshape(temp_coef,N,n_shift*n_theta);

        corr=temp_coef'* C(:,n_theta/2*k+1:end); % Compute the cross correlations with all other slices.
        corr=reshape(corr,n_shift*n_theta^2/2,n_proj-k);
        
        [temp_correlation, temp_idx] = max(real(corr)); % pick the largest one.
        ind1 = temp_correlation>correlation;
        
        correlation(ind1) = temp_correlation(ind1);
        idx(ind1) = temp_idx(ind1);
    end
    corrstack(k,k+1:n_proj)=correlation;
    [cl1, ~, cl2]=ind2sub([n_theta n_shift n_theta/2],idx);
    cl2(cl1>n_theta/2)=cl2(cl1>n_theta/2)+n_theta/2;
    cl1(cl1>n_theta/2)=cl1(cl1>n_theta/2)-n_theta/2;
    clstack(k,k+1:n_proj)=cl1;
    clstack(k+1:n_proj,k)=cl2';
end
