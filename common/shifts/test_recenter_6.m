% Test the effect of center of mass correction on the detection of common
% lines, for clean projections that are not centered (each projection is
% shifted randomly).
% The script:
%   1. Generate shifted clean projctions.
%   2. Finds common lines between the projections with and without taking
%      shifts into account.
%   3. Apply center of mass correction to all projections.
%   4. Verify that all common lines are correcctly detected, even without
%      searching over shifts (that is, all projections are correctly
%      aligned with respect to each other after applying center of mass
%      correction).
%
% Yoel Shkolnisky, November 2014.


clear;
K=10;     % Number of projections.
n=65;     % Size of each projection is nxn.
SNR=1e10; % Dummy SNR.
initstate;
[projs,~,~,rots]=cryo_gen_projections(n,K,SNR,3,1); % Generate shifted clean projections.
n_r=ceil(n/2);  % Radius of each projection
n_theta=360;    % Angular resolution of each projection (for computing its Fourier transform).
clstack_ref=clmatrix_cheat(rots,n_theta); % True common lines matrix.

pf=cryo_pft(projs,n_r,n_theta);
VERBOSE=0;

% Find common lines between projections while taking shifts into account.
open_log(0);
[clstack1,~,shift_equations1]=cryo_clmatrix(pf,K,VERBOSE,4,0.5);

% Percentage of correcctlydetected common lines.
prob1=comparecl(clstack1,clstack_ref,n_theta,5);
fprintf('Percentage of correct common lines when searching for shifts %7.4f\n',prob1);

% The last column of shift_equations holds the 1d relative shift of each
% pair of projections. Since the projectionsed are not centered, these
% relative shifts are not zero.
if nnz(shift_equations1(:,2*K+1))==0
    fprintf('Ref test OK - no 1d shifts between projections\n');
else
    figure;hist(shift_equations1(:,2*K+1));
    title('Relative 1d shifts between projections');
end

% Now search for common lines but ignore shifts between the projections.
[clstack2,~,shift_equations2]=cryo_clmatrix(pf,K,VERBOSE,0,1);
prob2=comparecl(clstack2,clstack_ref,n_theta,5);
fprintf('Percentage of correct common lines when ignoring shifts %7.4f\n',prob2);


% center all projections according to their center of mass.
projs_aligned=zeros(size(projs));
cmvec=zeros(K,2);
for j=1:K
    [projs_aligned(:,:,j),cm]=recenter(projs(:,:,j));
    cmvec(j,:)=cm;
end

% Find common lines among the corrected projections. First, while searching
% for shifts.
pf=cryo_pft(projs_aligned,n_r,n_theta);
open_log(0);
[clstack3,~,shift_equations3]=cryo_clmatrix(pf,K,VERBOSE,10,1);
prob3=comparecl(clstack3,clstack_ref,n_theta,5);
fprintf('Percentage of correct common lines after centering when looking for shifts %7.4f\n',prob3);

% Now find common lines between the corrected projections while not
% searching for shifts.
[clstack4,~,shift_equations4]=cryo_clmatrix(pf,K,VERBOSE,0,1);
prob4=comparecl(clstack4,clstack_ref,n_theta,5);
fprintf('Percentage of correct common lines after centering ignoring shifts %7.4f\n',prob4);

% As before, the last column of shift_equations3 indicates which pairs of
% projections have 1d shift along their common lines (meaning they are not
% centered relative to each other). There should be not such pairs, or very
% few, or the estimated relative shift is very small.
nz=nnz(shift_equations3(:,2*K+1));
if nz==0
    fprintf('Alignment fixed all shifts\n');
else
    fprintf('%d images pairs are not correctly aligned\n',nz);
    hist(shift_equations(:,2*K+1));
end
