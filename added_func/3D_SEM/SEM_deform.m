%% new method
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
%%
vmax=linspace(0,0.15,6);
n_theta=180;
shift_step = 1;
max_step = 15;
mask_radius = round(0.5*size(projections,1));
n_r = round(0.5*size(projections,1));
for i=1:6
load(['deform_vol_projection_' num2str(vmax(i)) '.mat'])
SEM_worker(data.projections,'SEM_deform_results.mat',1,...
    n_theta,n_r,max_step,shift_step,max_scale,data);
end