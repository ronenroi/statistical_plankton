function OUTPUT=fastrotate3z(INPUT,phi,M)
%FASTROTATE3Z Rotate a 3D volume around the z-axis.
% Input parameters:
%  INPUT    Volume to rotate, can be odd or even. 
%  phi      Rotation angle in degrees CCW. 
%  M        (Optional) Precomputed interpolation tables, as generated by
%           fastrotateprecomp. If M is given than phi is ignored. 
%
% Output parameters:
%  OUTPUT   The rotated volume.
%
% Examples:
%
%   rvol=fastrotate3z(vol,20);
%
%   M=fastrotateprecomp(size(vol,1),size(vol,2),20);
%   rvol=fastrotate(vol,[],M);
%
%Yoel Shkolnisky, November 2013.

[SzX, SzY, SzZ] =size(INPUT);

comp=0; % Determine if we need to compute the interpolation tables
if ~exist('M','var')
    comp=1;
elseif ~isstruct(M)
    comp=1;
elseif ~isfield(M,'Mx') || ~isfield(M,'My');
    comp=1;
end

if comp
    % Precompte M
    M=fastrotateprecomp(SzX,SzY,-phi);
end

OUTPUT=zeros(size(INPUT));
for k=1:SzZ
    im=squeeze(INPUT(:,:,k));
    rim=fastrotate(im,[],M);
    OUTPUT(:,:,k)=rim;
end