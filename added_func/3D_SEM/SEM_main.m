%% new method
clear;
close all
myPath = fileparts(mfilename('fullpath'));
cd ('../../');
initpath;
cd(myPath);
%%
vol_type=6;%0 simulation, 1 star, 2 Pyramimonas_longicauda 
switch vol_type
case 0
load('simulation_est/simulation_estimation_LUD_scale1.mat', 'data')
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = data.max_shift;
case 1 
dataFolder = '../../../data sets/jules images/daphnia'; projectionType = 'DarkField';  
projections = func_loadPlanktonData(dataFolder, projectionType);
max_shift = 15;
max_scale = 2;
data = -1;
case 2
dataFolder = '../../../data sets/Pyramimonas_longicauda'; projectionType = 'BrightField';  
projections = func_loadPlanktonData(dataFolder, projectionType);
max_shift = 15;
max_scale = 2;
data = -1;
case 3
dataFolder = '../../../data sets/mix4'; projectionType = 'BrightField';  
projections = func_loadPlanktonData_mix(dataFolder, projectionType);
max_shift = 15;
max_scale = 2;
data = -1;
case 4
dataFolder = '../../../data sets/oithona'; projectionType = 'DarkField';  
projections = func_loadPlanktonData(dataFolder, projectionType);
max_shift = 15;
max_scale = 2;
data = -1;
case 5
load('deform_vol_projection_0.06_3.mat', 'data')
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = 15;
case 6
%data = load(['simulations/simulated_data/scale' num2str(max_scale) '_new_rotations.mat']);
max_scale = data.max_scale;
projections = data.projections;
max_shift = 15;
end

%%
n_theta=180;
shift_step = 1;
mask_radius = round(0.75*size(projections,1));
n_r = round(0.75*size(projections,1));
SEM_worker_new(projections,'SEM_def0.06_3.mat',1,...
    n_theta,n_r,max_shift,shift_step,max_scale,data);