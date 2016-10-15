%% Reset environment
clear all;
close all;
clc;


inputFile = 'sphere2500vertigo.g2o';

%% Read in .g2o file
tic();
disp(sprintf('Loading file: %s ...', inputFile));
[measurements, edges_id, odomPoses] = readG2oDataset3D(inputFile, false);  
t = toc();
disp(sprintf('Processed input file in %g seconds\n', t));

%% Construct problem data structure
disp(sprintf('Constructing problem data matrices... '));
problem_data = construct_problem_data(measurements);

%% Set up optimization problem
r = 5;  % Maximum-rank parameter for Riemannian optimization
manopt_options.tolgradnorm = 1e-2;  % Stopping tolerance for Riemannian gradient
manopt_options.maxiter = 100;  % Maximum number of outer iterations (i.e. accepted update steps) to perform
manopt_options.maxinner = 300;  % Maximum number of iterations for the conjugate-gradient method used to compute approximate Newton steps
manopt_options.maxtime = 60*60;  % Maximum computation time to allow, in seconds

%% Solve Riemannian optimization problem
[Yopt, SDPLRval, info] =  Riemannian_optimization(problem_data, manopt_options, r);

%% Round solution and construct the estimate for the solution to the synchronization problem

disp(sprintf('\nRounding solution...\n'));
[Rhat, singular_values, determinants] = round_solution(Yopt, problem_data);
disp('Singular values of lifted solution Yopt:');
disp(singular_values);

% Compute primal rounded objective function value
RhatVal = evaluate_objective(Rhat, problem_data);

% Recover translational estimates
disp('Recovering translational estimates...');
t_hat = recover_translations(Rhat, problem_data);

%% Plot estimated poses
plot_loop_closures = false;
plot_poses_3D(t_hat, Rhat, edges_id, plot_loop_closures);
view(45, 25);

%% Output results
disp(sprintf('\nFINAL RESULTS:\n'));
disp(sprintf('Lifted solution value F(Yopt): %g', SDPLRval));
disp(sprintf('Rounded solution value F(Rhat): %g', RhatVal));
disp(sprintf('Gap: %g', RhatVal - SDPLRval));
disp(sprintf('Norm of Riemannian gradient at lifted solution (Yopt): %g', info.gradnorms(end)));
disp(sprintf('Total optimization time: %f seconds', info.total_elapsed_computation_time));

%% Compute minimum eigenvalue of Q - Lambda and corresponding eigenvector

% Compute the Lagrange multiplier corresponding to the projected solution
Lambda = compute_Lambda(Yopt, problem_data);



