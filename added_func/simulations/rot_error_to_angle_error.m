[MSE, err, O, aligned_rots1]=check_angle_err(aligned_rots_GSPR,data.ref_q);

% [MSE, err, O, aligned_rots1]=check_angle_err(aligned_rots,data.ref_q);
    

function [MSE, err, O, aligned_rots]=check_angle_err(est_inv_rots,q)
% The function checks the solution of SDP method in 
% http://www.math.princeton.edu/~amits/publications/common_lines_final_version.pdf
% The description of the SDP problem is (4.8) in the paper
% The result of the SDPLR is table 5.1 and table 5.3 in the paper
% 
% Input: 
%   est_inv_rots: the estimations of the orientations of the
%   images, a 3 \times 3 \times K matrix, where K is the number of images
%   and est_inv_rot_matrices(:,:,k) is the orientation of the k th matrix.
%
%   q: the true orientations of the images (quaternions), a 4 \times K matrix
%
% Output:
%   MSE: the error measurement
%   err: the errors of the estimated rotations.
%
% Based on Yoel and Amit's commonlines.m
% Lanhui Wang, Princeton University, Apr 6, 2011

K=size(q,2); % number of images
err=zeros(K,1);
est_rots_2=zeros(3,3,K);

%% Generate true orientations of images from q (return 'ref_inv_rot_matrices')
% Using the Euler formula
rot_matrices = zeros(3,3,K);

for k=1:K;
    rot_matrix = q_to_rot(q(:,k));
    rot_matrices(:,:,k) = rot_matrix;
end;

% calculate inverse rotation matrices (just transpose)
ref_inv_rot_matrices = zeros(3,3,K);


for k=1:K;
    rot_matrix = rot_matrices(:,:,k);
    inv_rot_matrix = rot_matrix'; % inv(R)=R^T
    ref_inv_rot_matrices(:,:,k) = inv_rot_matrix;
end;


corr_check1 = zeros(3);
corr_check2 = zeros(3);

for k=1:K;
    est_inv_rot = est_inv_rots(:,:,k);
    inv_rot = ref_inv_rot_matrices(:,:,k);
    corr_check1 = corr_check1 + inv_rot*est_inv_rot';
    est_inv_rot(:,3)=-est_inv_rot(:,3);%Try R3 = -R3
    est_inv_rot(3,:)=-est_inv_rot(3,:);
    est_rots_2(:,:,k)=est_inv_rot;
    corr_check2 = corr_check2 + inv_rot*est_inv_rot';
end;

MSE1 = 6 - 2 * sum(svd(corr_check1/K));
MSE2 = 6 - 2 * sum(svd(corr_check2/K));


if MSE1>MSE2
    MSE=MSE2;
    [S,~,D]=svd(corr_check2/K);
    est_rots=est_rots_2;
else
    MSE=MSE1;
    [S,~,D]=svd(corr_check1/K);
    est_rots=est_inv_rots;
end

O=S*D';

for k=1:K
%     err(k)=norm(ref_inv_rot_matrices(:,:,k)-O*est_rots(:,:,k),'fro');
gt(:,k) = [0 0 1] * ref_inv_rot_matrices(:,:,k);
est(:,k) = [0 0 1] * O*est_rots(:,:,k);
CosTheta = dot(gt(:,k),est(:,k))/(norm(gt(:,k))*norm(est(:,k)));
    ThetaInDegrees = real(acosd(CosTheta));
    if ThetaInDegrees>180
        ThetaInDegrees = 360-ThetaInDegrees;
    end
err(k)=ThetaInDegrees;
    aligned_rots(:,:,k) = O*est_rots(:,:,k);
end
plot_unit_sphere(gt,est)

end


function plot_unit_sphere(gt,est)
% Generate a unit sphere
theta=linspace(0,2*pi,20);
phi=linspace(0,pi,20);
[theta,phi]=meshgrid(theta,phi);
rho=1;
x=rho*sin(phi).*cos(theta);
y=rho*sin(phi).*sin(theta);
z=rho*cos(phi);
p1 = mesh(x,y,z);

hold on;
set(p1,'FaceAlpha',0);
p2 = plot3(gt(1,:),gt(2,:),gt(3,:),'ro');
set(p2,'MarkerFaceColor','red','LineWidth',2);
p3 = plot3(est(1,:),est(2,:),est(3,:),'go');
set(p2,'MarkerFaceColor','green','LineWidth',2);
end