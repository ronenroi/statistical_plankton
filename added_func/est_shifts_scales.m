function [ est_shifts, est_scales] = est_shifts_scales( shift_equations, scale_equations)
% Test the accuracy of the shift_equations/estimated shifts.
% 
% Input: 
%       shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%
%       ref_shifts: true shifts.
%
% Output:
%       est_shifts: estimated shifts.
%       err: l2 norm errors of the estimated shifts.
%
% Lanhui Wang, July 2, 2013
%    
K = (size(shift_equations, 2)-1)/2; % number of images
est_scales = [];

if exist('scale_equations')&&~isempty(scale_equations)
    % Magnification estimation using LS, and a l2 regularization on sizes of magnifications.
    % Regularization minimizes the energy of ||log(M)||^2
    est_scales=[scale_equations(:,1:end-1);eye(K)]...
        \[scale_equations(:,end);zeros(K,1)];
    est_scales = exp(est_scales);
end

% Shift estimation using LS, and a l2 regularization on sizes of shifts.
est_shifts=[shift_equations(:,1:end-1);eye(2*K)]...
    \[shift_equations(:,end);zeros(2*K,1)];
est_shifts = reshape(est_shifts, 2, K)';


end

