% This file contains a minimum working example demonstrating the use of the
% MATLAB distribution of SE-Sync, a certifiably correct algorithm for 
% synchronization over the special Euclidean group
%
% Copyright (C) 2016 by David M. Rosen

%% Reset environment
clear all;
close all;
clc;

%% Import SE-Sync
run('../import_SE_Sync.m');  % It's that easy :-)!

% Add helper code for loading files, plotting output, etc.
addpath(strcat(pwd, '/lib'));


%% Select dataset to run
% 3D datasets
prefix_3d = '../../data/3D/';

sphere2500 = strcat(prefix_3d, 'sphere2500vertigo');
sphere_a = strcat(prefix_3d,  'sphere_bignoise_vertex3');
torus = strcat(prefix_3d, 'torus3D');
grid = strcat(prefix_3d, 'grid3D');
garage = strcat(prefix_3d, 'parking-garage');
cubicle = strcat(prefix_3d, 'cubicle');

% 2D datasets
prefix_2d = '../../data/2D/';

m3500 =   strcat(prefix_2d, 'm3500');
city10000 = strcat(prefix_2d, 'city10000_vertigo');
intel = strcat(prefix_2d, 'intel_suger');
ais = strcat(prefix_2d, 'ais2klinik_suger');
eth = strcat(prefix_2d, 'ETHCampus_suger');

% Pick the dataset to run here
g2o_file = strcat(cubicle, '.g2o');

%% Read in .g2o file
tic();
disp(sprintf('Loading file: %s ...', g2o_file));
measurements = load_g2o_data(g2o_file);  
t = toc();
num_poses = max(max(measurements.edges));
num_measurements = length(measurements.kappa);
disp(sprintf('Processed input file %s in %g seconds', g2o_file, t));
disp(sprintf('Number of poses: %d', num_poses));
disp(sprintf('Number of measurements: %d\n', num_measurements));

%% Set Manopt options (if desired)
Manopt_opts.tolgradnorm = 1e-2;  % Stopping tolerance for norm of Riemannian gradient
Manopt_opts.rel_func_tol = 1e-5;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value
Manopt_opts.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxiter = 300;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
Manopt_opts.maxinner = 500;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
%manopt_options.maxtime = 60*60;  % Maximum computation time to allow, in seconds
%manopt_options.solver = @steepestdescent;  % Select Manopt solver to use: {trustregions (default), conjugategradient, steepestdescent}


%% Set SE-Sync options (if desired)
SE_Sync_opts.r0 = 5;  % Initial maximum-rank parameter at which to start the Riemannian Staircase
SE_Sync_opts.rmax = 10;  % Maximum maximum-rank parameter at which to terminate the Riemannian Staircase
SE_Sync_opts.eig_comp_rel_tol = 1e-4;  % Relative tolerance for the minimum-eigenvalue computation used to test for second-order optimality with MATLAB's eigs() function
SE_Sync_opts.min_eig_lower_bound = -1e-3;  % Minimum eigenvalue threshold for accepting a maxtrix as numerically positive-semidefinite

%% Run SE-Sync

% Pass explict settings for SE-Sync and/or Manopt
%[SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = SE_Sync(measurements, Manopt_opts, SE_Sync_opts);

% ... or ...

% Use default settings for everything
[SDPval, Yopt, xhat, Fxhat, se_sync_info, problem_data] = SE_Sync(measurements);

%% Plot resulting solution
plot_loop_closures = true;

if plot_loop_closures
    plot_poses(xhat.t, xhat.R, measurements.edges, '-b', .25);
else
    plot_poses(xhat.t, xhat.R);
end
axis tight;

%view(-90, 90);  % For plotting 3D but nearly-planar datasets

