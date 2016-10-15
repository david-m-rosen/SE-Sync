%% Reset environment
clear all;
close all;
clc;

%% Import SE-Sync toolbox
run('../import_SESync.m');  


%% Select dataset to run
% 3D datasets
is2D = false;
test = 'sphere2500vertigo';
%test ='sphere_bignoise_vertex3';
%test = 'torus3D';
%test = 'grid3D';
%test = 'parking-garage';
%test = 'cubicle';
%test = 'rim';

% 2D datasets
is2D = true;
%test = 'm3500_g2o';
%test = 'city10000_vertigo';
%test = 'intel_suger';
%test = 'ais2klinik_suger';
%test = 'ETHCampus_suger';
%test = 'kitti_00_odom_GT_LC.';
%test = 'kitti_02_odom_GT_LC.';
%test = 'kitti_05_odom_GT_LC';
%test = 'kitti_07_odom_GT_LC';
%test = 'kitti_08_odom_GT_LC';
%test = 'kitti_09_odom_GT_LC';

g2o_file = strcat(test, '.g2o');
workspace_file = strcat(test, '.mat'); 


%% Read in .g2o file
tic();
disp(sprintf('Loading file: %s ...', g2o_file));
measurements = construct_measurements(g2o_file, is2D);  
t = toc();
num_poses = max(max(measurements.edges));
num_measurements = length(measurements.kappa);
disp(sprintf('Processed input file %s in %g seconds', g2o_file, t));
disp(sprintf('Number of poses: %d', num_poses));
disp(sprintf('Number of measurements: %d\n', num_measurements));

%% Set Manopt options
manopt_options.tolgradnorm = 1e-2;  % Stopping tolerance for norm of Riemannian gradient
manopt_options.miniter = 1;  % Minimum number of outer iterations (i.e. accepted update steps) to perform
manopt_options.maxiter = 300;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
manopt_options.maxinner = 500;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
%manopt_options.maxtime = 60*60;  % Maximum computation time to allow, in seconds
%manopt_options.solver = @steepestdescent;  % Select Manopt solver to use: {trustregions (default), conjugategradient, steepestdescent}

%% Set SE-Sync options
se_sync_opts.r0 = 5;  % Initial maximum-rank parameter at which to start the Riemannian Staircase
se_sync_opts.rmax = 10;  % Maximum maximum-rank parameter at which to terminate the Riemannian Staircase
se_sync_opts.eig_comp_rel_tol = 1e-5;  % Relative tolerance for the minimum-eigenvalue computation used to test for second-order optimality with MATLAB's eigs() function
se_sync_opts.min_eig_threshold = -1e-4;  % Minimum eigenvalue threshold for accepting a maxtrix as numerically positive-semidefinite
se_sync_opts.relative_func_decrease_tol = 1e-6;  % Additional stopping criterion for Manopt: stop if the relative function decrease between two successive accepted iterates is less than this value

%% Run SE-Sync
[SDPval, Yopt, xhat, Fxhat, se_sync_info, auxiliary_data] = SE_Sync(measurements, manopt_options, se_sync_opts);

%% Save workspace
save(workspace_file);


%% Plot resulting solution
plot_loop_closures = true;

if plot_loop_closures
    plot_poses(xhat.t, xhat.R, measurements.edges, '-b', .25);
else
    plot_poses(xhat.t, xhat.R);
end
axis tight;

%View(-90, 30);  % For 3D but nearly-planar datasets


% Get rid of MATLAB's annoying border around the plotted stuff
set(gca,'position',[0 0 1 1],'units','normalized'); 

fig_filename = strcat(test, '.eps');
print(fig_filename, '-depsc');


