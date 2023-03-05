function [clstack,corrstack,shift_equations,scale_equations]=...
    cryo_clmatrix_cpu_roi(pf,max_shift,n_shifts,max_scale,Nscales ,rotations,sigma)
%
%
%   Generate common-lines matrix for the Fourier stack pf.
%
% Input parameters:
%   pf       3D array where each image pf(:,:,k) corresponds to the Fourier
%            transform of projection k.
%   NK       For each projection find its common-lines with NK other
%            projections. If NK is less than the total number a projection,
%            a random subset of NK projections is used. Default: n_proj.
%   verbose  Bitmask of debugging level.Bits:
%           0   silent
%           1   One line progress message (not written to log) (Default)
%           2   Print detailed debug messages
%           4   Draw common-line debugging plots
%           8   Draw shift estimation debugging plots
%   max_shift       Maximal 1D shift (in pixels)  to search between
%       common-lines. Default: 15.
%   shift_step      Resolution of shift estimation in pixels. Note that
%        shift_step can be any positive real number. Default: 1.
%   map_filter_radius      If nonzero, the common line between a pair
%       images is detected not by the pair of lines with the highest
%       correlation, but rather the pair of lines that both them and their
%       sorroundings given the best match. The radius for comparison is
%       determined by the value of map_filter_radius (Default 0).
%   ref_clmatrix    True common-lines matrix (for debugging).
%   ref_shifts_2d   True 2D shifts between projections (for debugging).
%
% Returned variables:
%
%   clstack     Common lines matrix. (k1,k2) and (k2,k1) contain the index
%       of the common line of projections k1 and k2. (k1,k2)  contains the
%       index of the common line in the projection k1. (k2,k1) contains the
%       index of the common line in k2.
%   corrstack   The correlation of the common line between projections k1
%       and k2. Since corrstack is symmetric, it contain entries only above
%       the diagonal. corrstack(k1,k2) measures how ''common'' is the between
%       projections k1 and k2. Small value means high-similariry.
%   shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%   shift_equations_map   2D array of size n_proj by n_proj. Entry (k1,k2)
%       is the index of the equation (row number) in the array
%       "shift_equations" that corresponds to the common line between
%       projections k1 and k2. shift_map is non-zero only for k1<k2.
%   clstack_mask  If ref_clmatrix is given, this array discribes which
%       common-lines were identified correcly. It is of size
%       n_projXn_projs, where entry (k1,k2) is 1 if the common-line between
%       projections k1 and k2 was correctly identified, and 0 otherwise.
%       This matrix will be non-zero only if bit 2 of verbose it set.
%
% Future version comment:
% The function assumes that the common-line is the pair of lines
% with maximum correlation. When the noise level is high, the
% maximal correlation usually does not correspond to the true
% common-line. In such cases, we would like to take several
% candidates for the common-line. See cryo_clmatrix_v4 for how
% to handle multiple candidates for common-line.
%
% Revisions:
%   02/03/09  Filename changed from cryo_clmatrix_v6.m to cryo_clmatrix.m.
%   02/04/09  Cleaning the code.
%   02/19/09  Maximal allowed deviation between common lines to be
%             considered as match was changed from 5 to 10 degrees.
%   02/24/09  Normalize the array pf only once (improves speed).
%   03/05/09  cryo_clmatrix_v3 created from cryo_clmatrix_v2
%   03/05/09  Exploit the conjugate symmetry of the Fourier transform to
%             compute correlations of length rmax instead of 2*rmax-1. This
%             gives a factor of 2 in performance.
%   15/06/13  Replace repmat with bsxfnu when comparing shifted lines. This
%             should be faster and makes the code easier to port to GPU.
%   28/2/17   Rename cryo_clmatrix to cryo_clmatrix_cpu.

msg=[];

T=size(pf,2);

if mod(T,2)~=0
    error('n_theta must be even');
end

% pf is of size n_rxn_theta. Convert pf into an array of size
% (2xn_r-1)xn_theta, that is, take then entire ray through the origin, but
% thake the angles only up PI.
% This seems redundant: The original projections are real, and thus
% each ray is conjugate symmetric. We therefore gain nothing by taking
% longer correlations (of length 2*n_r-1 instead of n_r), as the two halfs
% are exactly the same. Taking shorter correlation would speed the
% computation by a factor of two.

[n_r,n_theta,n_proj]=size(pf);
shift_step = 2*max_shift/(n_shifts-1);
scales = exp(linspace(-log(max_scale), log(max_scale), Nscales));

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
pf3=reshape(pf,N,n_theta,n_proj);



clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.


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
for k=1:n_theta
    for l=1:n_proj
        pf3(:,k,l)=pf3(:,k,l)/norm(pf3(:,k,l));
    end
    
end




rk2=rk(1:rmax);
for k1=1:n_proj
    
    proj1=pf3(:,:,k1);
    % P1=proj1(1:rmax,:);  % Take half ray plus the DC
    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    %     if norm(proj1(rmax+1,:))>1.0e-13
    %         error('DC component of projection is not zero');
    %     end
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(proj1, size(proj1).*[M, 1]);
        if M<1.0
            scale_coeff = zeros(size(proj1));
            scale_coeff(end-size(scaled_mag_coeff,1)+1:end,:)=scaled_mag_coeff;
        else
            scale_coeff = scaled_mag_coeff(end-N+1:end,:);
        end
        
        for kk=1:size(scale_coeff,2)
            scale_coeff(:,kk) = (1/norm(scale_coeff(:,kk)))*scale_coeff(:,kk);
        end
        %P1_flipped=conj(P1);
        
        for k2=k1+1:n_proj
            Ri=rotations(:,:,k1);
            Rj=rotations(:,:,k2);
            [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
            
            cl_i_vec = mod(cij-sigma : cij+sigma,n_theta/2)+1 ;
            cl_j_vec = mod(cji-sigma : cji+sigma,n_theta)+1 ;
            P1=scale_coeff(:,cl_i_vec);  %
            proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
            P2=proj2(:,cl_j_vec);
            
            %         if norm(proj2(rmax+1,:))>1.0e-13
            %             error('DC component of projection is not zero');
            %         end
            
            % Find the shift that gives best correlation.
            for shiftidx=1:n_shifts
                shift=-max_shift+(shiftidx-1)*shift_step;
                shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1));
                %shift_phases=repmat(shift_phases,1,n_theta);
                
                % No need to renormalize proj1_shifted and
                % proj1_shifted_flipped since multiplication by phases
                % does not change the norm, and proj1 is already normalized.
                %P1_shifted=P1.*shift_phases;
                %P1_shifted_flipped=P1_flipped.*shift_phases;
                
                P1_shifted=bsxfun(@times,P1,shift_phases);
                %             P1_shifted_flipped=bsxfun(@times,P1_flipped,shift_phases);
                
                % Compute correlations in the positive r direction
                C=real(P1_shifted'*P2);
                
                % Compute correlations in the negative r direction
                %             C2=2*real(P1_shifted_flipped'*P2);
                
                
                %             C = [C1,C2];
                
                
                
                [sval,sidx]=max(C(:));
                
                % correlation than previously known.
                
                if sval>corrstack(k1,k2)
                    %                 [cl1,cl2]=ind2sub([n_theta 2*n_theta],sidx);
                    %                 clstack(k1,k2)=cl1;
                    %                 clstack(k2,k1)=cl2;
                    %                 corrstack(k1,k2)=sval;
                    %                 shifts_1d(k1,k2)=shift;
                    %                 improved_correlation=1;
                    [cli,clj]=ind2sub([length(cl_i_vec) length(cl_i_vec)],sidx);
                    cl1 = cl_i_vec(cli);
                    cl2 = cl_j_vec(clj);
                    cl2(cl1>n_theta/2)=cl2(cl1>n_theta/2)+n_theta/2;
                    cl1(cl1>n_theta/2)=cl1(cl1>n_theta/2)-n_theta/2;
                    clstack(k1,k2)=cl1;
                    clstack(k2,k1)=cl2;
                    corrstack(k1,k2)=sval;
                    shifts_1d(k1,k2)=shift;
                    scale_1d(k1,k2) = log(M);
                end
                
                
            end
            
            
        end
    end
    
    
end
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

scale_equations=sparse(scale_I(1:2*shift_equation_idx),...
    scale_J(1:2*shift_equation_idx),scale_eq(1:2*shift_equation_idx),...
    shift_equation_idx, n_proj);
scale_equations=[scale_equations scale_b(1:shift_equation_idx)];
end