function [clstack,corrstack,shift_equations,shift_equations_map,clstack_mask,scale_equations]=...
    cryo_clmatrix_cpu_roi_scale(pf,max_shift,n_shifts,max_scale,Nscales ,rotations,sigma)
%
%

msg=[];

T=size(pf,2);

if mod(T,2)~=0
    error('n_theta must be even');
end

%pf=[flipdim(pf(2:end,T/2+1:end,:),1) ; pf(:,1:T/2,:) ];

scales = exp(linspace(-log(max_scale), log(max_scale), Nscales));

[n_r, n_theta, n_proj]=size(pf);


%% Check input parameters and set debug flags.
NK=n_proj;
if (nargin<2) || (NK==-1)
    NK=n_proj; % Number of common-line pairs to compute for each projection
end

if ~exist('verbose','var')
    verbose=1;
end

if ~exist('max_shift','var')
    max_shift=15; % Maximal shift between common-lines in pixels. The
    % shift  is from -max_shift to max_shift.
end

if ~exist('shift_step','var')
    shift_step=1.0; % Resolution of shift estimation in pixels.
end
shift_step = 2*max_shift/(n_shifts-1);

if ~exist('map_filter_radius','var')
    map_filter_radius=0;
end

if ~exist('ref_clmatrix','var') || isempty(ref_clmatrix)
    ref_clmatrix=0;
end

if ~exist('ref_shifts_2d','var') || isempty(ref_shifts_2d)
    ref_shifts_2d=0;
end


% Set flag for progress and debug messages
verbose_progress=0;
verbose_detailed_debugging=0;
verbose_plot_cl=0;
verbose_plot_shifts=0;

if bitand(verbose,1)
    verbose_progress=1;
end

found_ref_clmatrix=0;
if ~isscalar(ref_clmatrix)
    found_ref_clmatrix=1;
else
    if verbose>0
        log_message('Reference clmatrix not found');
    end
end

found_ref_shifts=0;
if ~isscalar(ref_shifts_2d)
    found_ref_shifts=1;
else
    if verbose>0
        log_message('Reference shifts not found');
    end
end

if bitand(verbose,2)
    verbose_detailed_debugging=1;
    verbose_progress=0;
end

if bitand(verbose,4)
    if isscalar(ref_clmatrix)
        log_message('Common-lines plots not available. Reference clmatrix is missing\n');
    end
    verbose_plot_cl=1;
end

if bitand(verbose,8)
    if isscalar(ref_clmatrix) || isscalar(ref_shifts_2d)
        log_message('Only partial information will be plotted. Reference clmatrix or shifts are missing\n');
    end
    verbose_plot_shifts=1;
end

if verbose~=0
    log_message('Verbose mode=%d',verbose);
end

%%

clstack=zeros(n_proj,n_proj);      % Common lines-matrix.
corrstack=zeros(n_proj,n_proj);    % Correlation coefficient for each common-line.
clstack_mask=zeros(n_proj,n_proj); % Which common-lines were correctly identified.

refcorr=zeros(n_proj,n_proj); % Correlation between true common-lines.
thetadiff=zeros(n_proj,n_proj); % Angle between true and estimated common lines.

%% Allocate variables used for shift estimation

shifts_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines.

ref_shifts_1d=zeros(n_proj,n_proj); % True shift along the common-line
% between each pair of projections. Computed from the reference 2D
% shifts.

shift_estimation_error=zeros(n_proj,n_proj); % The difference between the
% estimated shift along each common line and the true shift.

% Based on the estimated common-lines, construct the equations for
% determining the 2D shift of each projection. The shift equations are
% represented using a sparse matrix, since each row in the system contains
% four non-zeros (as it involves exactly four unknowns).
% The variables below are used to construct this sparse system. The k'th
% non-zero element of the equations matrix is stored at index
% (shift_I(k),shift_J(k)).
shift_I=zeros(4*n_proj*NK,1);  % Row index for sparse equations system.

shift_J=zeros(4*n_proj*NK,1);  % Column index for sparse equations system.

shift_eq=zeros(4*n_proj*NK,1); % The coefficients of the center estimation
% system ordered as a single vector.

shift_equations_map=zeros(n_proj); % Entry (k1,k2) is the index of the
% euqation for the common line of projections k1 and k2.

shift_equation_idx=1;  % The equation number we are currently processing.
shift_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system.
dtheta=pi/n_theta; % Not 2*pi/n_theta, since we divided n_theta by 2 to
% take rays of length 2*n_r-1.
%% Allocate variables used for magnification estimation

scale_1d=zeros(n_proj,n_proj);     % Estimated 1D shift between common-lines.
scale_I=zeros(2*n_proj^2,1);  % Row index for sparse equations system.

scale_J=zeros(2*n_proj^2,1);  % Column index for sparse equations system.

scale_eq=zeros(2*n_proj^2,1); % The coefficients of the center estimation
% system ordered as a single vector.

scale_b=zeros(n_proj*(n_proj-1)/2,1);   % Right hand side of the system
if verbose>0
    log_message('Shift estimation parameters: max_shift=%d   shift_step=%d',max_shift,shift_step);
end

if verbose>0
    log_message('map_filter_radius = %d',map_filter_radius);
end


%% Debugging handles and variables

matched_cl=0;  % How many times the estimated common-line is close (to
% within a prescribed tolerance) to the true common-line.
% Used for debugging.



%% Search for common lines between pairs of projections

% Construct filter to apply to each Fourier ray.
rmax=(size(pf,1)-1);
rk=0:rmax; rk=rk(:);
H=sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2));
H=repmat(H(:),1,n_theta);  % Filter for common-line detection.

% Bandpass filter and normalize each ray of each projection.
% XXX We do not override pf since it is used to debugging plots below. Once
% XXX these debugging plots are removed, replace pf3 by pf. This will save
% XXX a lot of memory.
pf3=pf;
for k=1:n_proj
    proj=pf(:,:,k);
    proj=proj.*H;
    proj(1:2,:)=0;
    proj=cryo_raynormalize(proj);
    pf3(:,:,k)=proj;
end

rk2=rk(1:rmax);
for k1=1:n_proj
    
    proj1=pf3(:,:,k1);
    
    
    % Make sure the DC component is zero. This is assumed  below in
    % computing correlations.
    if norm(proj1(1,:))>1.0e-13
        error('DC component of projection is not zero');
    end
    for scale_index=1:Nscales
        M=scales(scale_index);
        scaled_mag_coeff = imresize(proj1, size(proj1).*[M, 1]);
        if M<1.0
            scale_coeff = zeros(size(proj1));
            scale_coeff(end-size(scaled_mag_coeff,1)+1:end,:)=scaled_mag_coeff;
        else
            scale_coeff = scaled_mag_coeff(end-n_r+1:end,:);
        end
        
        for kk=1:size(scale_coeff,2)
            scale_coeff(:,kk) = (1/norm(scale_coeff(:,kk)))*scale_coeff(:,kk);
        end
        for k2=k1+1:n_proj
            Ri=rotations(:,:,k1);
            Rj=rotations(:,:,k2);
            [cij,cji]=commonline_R(Ri.',Rj.',n_theta);
            
            cl_i_vec = mod(cij-sigma : cij+sigma,n_theta)+1 ;
            cl_j_vec = mod(cji-sigma : cji+sigma,n_theta)+1 ;
            P1=scale_coeff(1:rmax,cl_i_vec);  % Take half ray plus the DC
            %P1_flipped=conj(P1);
            proj2=pf3(:,:,k2); % proj1 and proj2 are both normalized to unit norm.
            P2=proj2(1:rmax,cl_j_vec);
            
            if norm(proj2(1,:))>1.0e-13
                error('DC component of projection is not zero');
            end
            
            
            
            
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
                %P1_shifted_flipped=bsxfun(@times,P1_flipped,shift_phases);
                
                % Compute correlations in the positive r direction
                C1=2*real(P1_shifted'*P2);
                
                % Compute correlations in the negative r direction
                %C2=2*real(P1_shifted_flipped'*P2);
                
                
                % C = [C1,C2];
                C=C1;
                if map_filter_radius > 0
                    C = cryo_average_clmap(C, map_filter_radius);
                end
                
                [sval,sidx]=max(C(:));
                
                improved_correlation=0; % Indicates that we found a better
                % correlation than previously known.
                
                if sval>corrstack(k1,k2)
                    %                 [cl1,cl2]=ind2sub([n_theta 2*n_theta],sidx);
                    [cl1,cl2]=ind2sub([length(cl_i_vec) length(cl_i_vec)],sidx);
                    
                    clstack(k1,k2)=cl_i_vec(cl1);
                    clstack(k2,k1)=cl_j_vec(cl2);
                    
                    %                 if cl2>length(cl_i_vec)
                    %                 clstack(k2,k1)=cl_j_vec(cl2-length(cl_i_vec))+n_theta;
                    %                 end
                    
                    corrstack(k1,k2)=sval;
                    shifts_1d(k1,k2)=shift;
                    scale_1d(k1,k2) = log(M);
                    improved_correlation=1;
                end
                
                
            end
            
            
            
            
            
            
            
            
        end
        %                         % Create a shift equation for the projections pair (k1,k2).
        %                 idx=4*(shift_equation_idx-1)+1:4*shift_equation_idx;
        %                 shift_alpha=(clstack(k1,k2)-1)*dtheta;  % Angle of common ray in projection 1.
        %                 shift_beta= (clstack(k2,k1)-1)*dtheta;  % Angle of common ray in projection 2.
        %                 shift_I(idx)=shift_equation_idx; % Row index to construct the sparse equations.
        %                 shift_J(idx)=[2*k1-1 2*k1 2*k2-1 2*k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        %                 shift_b(shift_equation_idx)=shifts_1d(k1,k2); % Right hand side of the current equation.
        %
        %                 % Compute the coefficients of the current equation.
        %                 if shift_beta<pi
        %                     shift_eq(idx)=[sin(shift_alpha) cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        %                 else
        %                     shift_beta=shift_beta-pi; % In the derivation we assume that all angles are less
        %                     % than PI where angles larger than PI are assigned
        %                     % nigative orientation.
        %                     shift_eq(idx)=[-sin(shift_alpha) -cos(shift_alpha) -sin(shift_beta) -cos(shift_beta)];
        %                 end
        %
        %                 shift_equations_map(k1,k2)=shift_equation_idx;  % For each pair (k1,k2), store the index of its equation.
        %
        %                 % Scale equations
        %                 sidx=2*(shift_equation_idx-1)+1:2*shift_equation_idx;
        %                 scale_I(sidx)=shift_equation_idx; % Row index to construct the sparse equations.
        %                 scale_J(sidx)=[k1 k2]; % Columns of the shift variables that correspond to the current pair (k1,k2).
        %                 scale_b(shift_equation_idx)=scale_1d(k1,k2); % Right hand side of the current equation.
        %                 scale_eq(sidx)=[1 -1];
        %                                 shift_equation_idx=shift_equation_idx+1;
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

% corrstack(corrstack~=0)=1-corrstack(corrstack~=0);
%
% % Construct least-squares for the two-dimensioal shifts.
% shift_equation_idx=shift_equation_idx-1;
% shift_equations=sparse(shift_I(1:4*shift_equation_idx),...
%     shift_J(1:4*shift_equation_idx),shift_eq(1:4*shift_equation_idx),...
%     shift_equation_idx,2*n_proj);
%
% shift_equations=[shift_equations shift_b(1:shift_equation_idx)];
%
% scale_equations=sparse(scale_I(1:2*shift_equation_idx),...
%     scale_J(1:2*shift_equation_idx),scale_eq(1:2*shift_equation_idx),...
%     shift_equation_idx, n_proj);
% scale_equations=[scale_equations scale_b(1:shift_equation_idx)];
end




