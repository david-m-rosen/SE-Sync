function [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = SE_Sync(measurements, Manopt_opts, SE_Sync_opts, Y0)
%function [SDPval, Yopt, xhat, Fxhat, SE_Sync_info, problem_data] = SE_Sync(measurements, Manopt_opts, SE_Sync_opts, Y0)
%
% SE-Sync: A certifiably correct algorithm for synchronization over the
% special Euclidean group
%
%
% INPUTs:
%
% measurements:  A MATLAB struct containing the data describing the special
%   Euclidean synchronization problem (see eq. (11) in the paper for
%   details). Specifically, measurements must contain the following fields:
%   edges:  An (mx2)-dimensional matrix encoding the edges in the measurement
%     network; edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform x_i^{-1} x_j.  NB:  This indexing scheme requires
%     that the states x_i are numbered sequentially as x_1, ... x_n.
%   R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
%   t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
%   kappa:  An m-dimensional cell array whose kth element gives the
%     precision of the rotational part of the kth measurement.
%   tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.
%
% Manopt_opts [optional]:  A MATLAB struct containing various options that
%       determine the behavior of Manopt's Riemannian truncated-Newton
%       trust-region method, which we use to solve instances of the
%       rank-restricted form of the semidefinite relaxation.  This struct
%       contains the following [optional] fields (among others, see the
%       Manopt documentation)
%   tolgradnorm:  Stopping criterion; norm tolerance for the Riemannian gradient
%   rel_func_tol:  An additional stopping criterion for the Manopt
%     solver.  Terminate whenever the relative decrease in function value
%     between subsequenct iterations is less than this value (in the range
%     (0,1) ).
%   maxinner:  Maximum number of Hessian-vector products to evaluate as part
%      of the truncated conjugate-gradient procedure used to compute update
%      steps.
%   miniter:  Minimum number of outer iterations (update steps).
%   maxiter:  Maximum number of outer iterations (update steps).
%   maxtime:  Maximum permissible elapsed computation time (in seconds).
%   preconditioner: A string specifying an (optional) preconditioner to use
%      when computing inexact Newton steps via the truncated conjugate
%      gradient algorithm.  Possible values are:
%      - 'ichol':  Use a zero-fill incomplete Cholesky preconditioner [default]
%      - 'Jacobi':  Use a simple diagonal Jacobi preconditioner
%      - 'none':  Do not precondition  
%
% SE_Sync_opts [optional]:  A MATLAB struct determining the behavior of the
%       SE-Sync algorithm.  This struct contains the following [optional]
%       fields:
%   r0:  The initial value of the maximum-rank parameter r at which to
%      start the Riemannian Staircase
%   rmax:  The maximum value of the maximum-rank parameter r.
%   eig_comp_max_iters:  The maximum number of Lanczos iterations to
%      perform when computing the minimum eigenvalue
%   min_eig_num_tol:  Lower bound for the minimum eigenvalue in 
%      order to consider the matrix Q - Lambda to be positive semidefinite.
%      Typical values here should be small-magnitude numbers, e.g. 10^-4
%   Cholesky:  A Boolean value indicating whether to compute orthogonal
%      projections onto the cycle space of G using a cached Cholesky
%      factorization of Ared*Ared' or by applying an orthogonal (QR)
%      decomposition.  The former method may be faster on smaller problems, 
%      but the latter is more numerically stable [default: false]
%   init:  A string specifying the initialization procedure to use if no
%      initial point Y0 is passed.  Options are 'chordal' or 'random'.  If
%      no option is specified, 'chordal' is used as a default
%
% Y0:  [Optional]  An initial point on the manifold St(d, r)^n at which to
%      initialize the first Riemannian optimization problem.  If this
%      parameter is not passed, a randomly-sampled point is used instead.
%
%
% OUTPUTS:
%
% SDPval:  The optimal value of the semidefinite relaxation
% Yopt:  A symmetric factor of an optimal solution Zopt = Yopt' * Yopt for
%      the semidefinite relaxation.
% xhat: A struct containing the estimate for the special Euclidean
%   synchronization problem.  It has the following two fields:
%   Rhat:  A d x dn matrix whose (dxd)-block elements give the rotational
%   state estimates.
%   that: a d x n matrix whose columsn give the translational state estimates.
% Fxhat:  The objective value of the rounded solution xhat.
%
% SE_Sync_info:  A MATLAB struct containing various possibly-interesting
%   bits of information about the execution of the SE-Sync algorithm.  The
%   fields are:
%   mat_contruct_times:  The elapsed time needed to construct the auxiliary
%     system matrices contained in 'problem_data'
%   init_time:  The elapsed time needed to compute the initial point for
%     the Riemannian Staircase.
%   optimization_times:  A vector containing the elapsed computation times
%     for solving the optimization problem at each level of the Riemannian
%     Staircase.
%   SDPLRvals:  A vector containing the optimal value of the optimization
%     problem solved at each level of the Riemannian Staircase
%   min_eig_times:  A vector containing the elapsed computation times for
%     performing the minimum-eigenvalue computation necessary to check for
%     optimality of Yopt as a solution of the SDP after solving each
%     Riemannian optimization problem to first-order.
%   min_eig_vals:  A vector containing the corresponding minimum
%      eigenvalues.
%   Yvals:  A cell array containing the sequence of iterates obtained by
%      the Riemannian Staircase
%   gradnorms:  Norms of the gradients at the sequence of iterates obtained
%      by the Riemannian Staircase
%   total_computation_time:  The elapsed computation time of the complete
%      SE-Sync algorithm
%   manopt_info:  The info struct returned by the Manopt solver for the
%      during its last execution (i.e. when solving the last explored level
%      the Riemannian Staircase).
%
% problem_data:  A MATLAB struct containing several auxiliary matrices
% constructed from the input measurements that are used internally
% throughout the SE-Sync algorithm.  Specifically, this struct contains the
% following fields:
%   n:  The number of group elements (poses) to estimate
%   m:  The number of relative measurements
%   d:  The dimension of the Euclidean space on which these group elements
%       act (generally d is 2 or 3).
%   LWtau:  The Laplacian for the translational weight graph W^tau.
%   ConLap:  The connection Laplacian for the set of rotational
%       measurements; see eq. (15) in the paper.
%   A:  An oriented incidence matrix for the directed graph of
%       measurements; see eq. (7) in the paper
%   Ared:  The reduced oriented incidence matrix obtained by removing the
%       final row of A.
%   L:  A sparse lower-triangular Cholesky factor of the reduced Laplacian
%       of the translational weight graph
%   T:  The sparse matrix of translational observations defined in eq. (24)
%       in the paper
%   Omega:  The diagonal matrix of translational measurement precisions;
%       defined in eq. (23).
%   V:  The sparse translational data matrix defined in eq. (16) in the
%       paper.

% Copyright (C) 2016, 2017 by David M. Rosen


fprintf('\n\n========== SE-Sync ==========\n\n');

timerVal = tic();


%% INPUT PARSING

% SE-Sync settings:
fprintf('ALGORITHM SETTINGS:\n\n');

if nargin < 3
    disp('Using default settings for SE-Sync:');
    SE_Sync_opts = struct;  % Create empty structure
else
    disp('SE-Sync settings:');
end

if isfield(SE_Sync_opts, 'r0')
    fprintf(' Initial level of Riemannian Staircase: %d\n', SE_Sync_opts.r0);
else
    SE_Sync_opts.r0 = 5;
    fprintf(' Setting initial level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.r0);
end

if isfield(SE_Sync_opts, 'rmax')
    fprintf(' Final level of Riemannian Staircase: %d\n', SE_Sync_opts.rmax);
else
    SE_Sync_opts.rmax = 7;
    fprintf(' Setting final level of Riemannian Staircase to %d [default]\n', SE_Sync_opts.rmax);
end

% if isfield(SE_Sync_opts, 'eig_comp_rel_tol')
%     fprintf(' Relative tolerance for minimum eigenvalue computation in test for positive semidefiniteness: %g\n', SE_Sync_opts.eig_comp_rel_tol);
% else
%     SE_Sync_opts.eig_comp_rel_tol = 1e-8;
%     fprintf(' Setting relative tolerance for minimum eigenvalue computation in test for positive semidefiniteness to: %g [default]\n', SE_Sync_opts.eig_comp_rel_tol);
% end

if isfield(SE_Sync_opts, 'eig_comp_max_iters')
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g\n', SE_Sync_opts.eig_comp_max_iters);
else
    SE_Sync_opts.eig_comp_max_iters = 2000;
    fprintf(' Maximum number of iterations to perform for minimum eigenvalue computation in test for positive semidefiniteness: %g [default]\n', SE_Sync_opts.eig_comp_max_iters);
end

if isfield(SE_Sync_opts, 'min_eig_num_tol')
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g\n', SE_Sync_opts.min_eig_num_tol);
else
    SE_Sync_opts.min_eig_num_tol = 1e-4;
    fprintf(' Tolerance for accepting an eigenvalue as numerically nonnegative in optimality verification: %g [default]\n', SE_Sync_opts.min_eig_num_tol);
end


if ~isfield(SE_Sync_opts, 'Cholesky')
    fprintf(' Using QR decomposition to compute orthogonal projection [default]\n');
    SE_Sync_opts.Cholesky = false;
else
    if SE_Sync_opts.Cholesky
        fprintf(' Using Cholesky decomposition to compute orthogonal projection\n');
    else
        fprintf(' Using QR decomposition to compute orthogonal projection\n');
    end
end

if ~isfield(SE_Sync_opts, 'init')
    fprintf(' Initialization method: chordal [default]\n');
    SE_Sync_opts.init = 'chordal';
else
    if strcmp(SE_Sync_opts.init, 'chordal')
        fprintf(' Initialization method: chordal\n');
    elseif strcmp(SE_Sync_opts.init, 'random');
        fprintf(' Initialization method: random\n');
    else
        error(sprintf('Initialization option "%s" not recognized!  (Supported options are "chordal" or "random"\n', SE_Sync_opts.init));
    end
end



fprintf('\n');

%% Manopt settings:

if nargin < 2
    disp('Using default settings for Manopt:');
    Manopt_opts = struct;  % Create empty structure
else
    disp('Manopt settings:');
end

if isfield(Manopt_opts, 'tolgradnorm')
    fprintf(' Stopping tolerance for norm of Riemannian gradient: %g\n', Manopt_opts.tolgradnorm);
else
    Manopt_opts.tolgradnorm = 1e-2;
    fprintf(' Setting stopping tolerance for norm of Riemannian gradient to: %g [default]\n', Manopt_opts.tolgradnorm);
end

if isfield(Manopt_opts, 'rel_func_tol')
    fprintf(' Stopping tolerance for relative function decrease: %g\n', Manopt_opts.rel_func_tol);
else
    Manopt_opts.rel_func_tol = 1e-5;
    fprintf(' Setting stopping tolerance for relative function decrease to: %g [default]\n', Manopt_opts.rel_func_tol);
end

if isfield(Manopt_opts, 'maxinner')
    fprintf(' Maximum number of Hessian-vector products to evaluate in each truncated Newton iteration: %d\n', Manopt_opts.maxinner);
else
    Manopt_opts.maxinner = 1000;
    fprintf(' Setting maximum number of Hessian-vector products to evaluate in each truncated Newton iteration to: %d [default]\n', Manopt_opts.maxinner);
end

if isfield(Manopt_opts, 'miniter')
    fprintf(' Minimum number of trust-region iterations: %d\n', Manopt_opts.miniter);
else
    Manopt_opts.miniter = 1;
    fprintf(' Setting minimum number of trust-region iterations to: %d [default]\n', Manopt_opts.miniter);
end

if isfield(Manopt_opts, 'maxiter')
    fprintf(' Maximum number of trust-region iterations: %d\n', Manopt_opts.maxiter);
else
    Manopt_opts.maxiter = 500;
    fprintf(' Setting maximum number of trust-region iterations to: %d [default]\n', Manopt_opts.maxiter);
end

if isfield(Manopt_opts, 'maxtime')
    fprintf(' Maximum permissible elapsed computation time [sec]: %g\n', Manopt_opts.maxtime);
end

if ~isfield(Manopt_opts, 'preconditioner')
    fprintf(' Using incomplete zero-fill Cholesky preconditioner for truncated conjugate gradient inexact Newton step computations [default]\n');
    Manopt_opts.preconditioner = 'ichol';
else
    if(strcmp(Manopt_opts.preconditioner, 'ichol'))
        fprintf(' Using incomplete zero-fill Cholesky preconditioner for truncated conjugate gradient inexact Newton step computations\n');
    elseif(strcmp(Manopt_opts.preconditioner, 'Jacobi'))
        fprintf(' Using Jacobi preconditioner for truncated conjugate gradient inexact Newton step computations\n');
    elseif(strcmp(Manopt_opts.preconditioner, 'none'))
        fprintf(' Using unpreconditioned truncated conjugate gradient for inexact Newton step computations\n');
    else
        error(sprintf('Initialization option "%s" not recognized!  (Supported options are "Jacobi" or "none"\n', Manopt_opts.preconditioner));
    end
end

        






%% Construct problem data matrices from input
fprintf('\n\nINITIALIZATION:\n\n');
disp('Constructing auxiliary data matrices from raw measurements...');
aux_time_start = tic();
problem_data = construct_problem_data(measurements);
auxiliary_matrix_construction_time = toc(aux_time_start);
fprintf('Auxiliary data matrix construction finished.  Elapsed computation time: %g seconds\n\n', auxiliary_matrix_construction_time);

%% Construct (Euclidean) preconditioning function handle, if desired
if isfield(Manopt_opts, 'preconditioner')
    precon_construction_start_time = tic();
    
    if (strcmp(Manopt_opts.preconditioner, 'ichol'))
        fprintf('Constructing incomplete Cholesky preconditioner... ');
        
        LGrho = problem_data.ConLap;
        
        % Regularize this matrix by adding a very small positive
        % multiple of the identity to account for the fact that the
        % rotational connection Laplacian is singular
        ichol_opts.diagcomp = 1e-3;
        
        % Compute incomplete zero-fill
        L = ichol(LGrho, ichol_opts);
        LT = L';
        
        precon = @(u) LT \ (L \ u);
        
    elseif(strcmp(Manopt_opts.preconditioner, 'Jacobi'))
        fprintf('Constructing Jacobi preconditioner... ');
        J = problem_data.ConLap;
        
        % Extract diagonal elements
        D = spdiags(J, 0);
        
        % Invert these
        Dinv = 1 ./ D;
        
        % Construct diagonal matrix with this size
        Pinv = spdiags(Dinv, 0, problem_data.d * problem_data.n, problem_data.d * problem_data.n);
        
        % Set preconditioning function
        precon = @(u) Pinv * u;
    end
    
    if(~strcmp(Manopt_opts.preconditioner, 'none'))
        precon_construction_end_time = toc(precon_construction_start_time);
        fprintf('Elapsed computation time: %g seconds\n\n', precon_construction_end_time);
    end
end



%% INITIALIZATION

% The maximum number of levels in the Riemannian Staircase that we will
% need to explore
max_num_iters = SE_Sync_opts.rmax - SE_Sync_opts.r0 + 1;

% Allocate storage for state traces
optimization_times = zeros(1, max_num_iters);
SDPLRvals = zeros(1, max_num_iters);
min_eig_times = zeros(1, max_num_iters);
min_eig_vals = zeros(1, max_num_iters);
gradnorms = [];
Yvals = {};

% Set up Manopt problem

% We optimize over the manifold M := St(d, r)^N, the N-fold product of the
% (Stiefel) manifold of orthonormal d-frames in R^r.
manopt_data.M = stiefelstackedfactory(problem_data.n, problem_data.d, SE_Sync_opts.r0);

% Check if an initial point was supplied
if nargin < 4
    if strcmp(SE_Sync_opts.init, 'chordal')
        fprintf('Computing chordal initialization...\n');
        init_time_start = tic();
        Rchordal = chordal_initialization(measurements);
        Y0 = vertcat(Rchordal, zeros(SE_Sync_opts.r0 - problem_data.d, problem_data.d * problem_data.n));
        init_time = toc(init_time_start);
    else  % Use randomly-sampled initialization
        fprintf('Randomly sampling an initial point on St(%d,%d)^%d ...\n', problem_data.d, SE_Sync_opts.r0, problem_data.n);
        init_time_start = tic();
        % Sample a random point on the Stiefel manifold as an initial guess
        Y0 = manopt_data.M.rand()';
        init_time = toc(init_time_start);
    end
    fprintf('Elapsed computation time: %g seconds\n', init_time);
else
    fprintf('Using user-supplied initial point Y0 in Riemannian Staircase\n\n');
    init_time = 0;
end

% Check if a solver was explicitly supplied
if(~isfield(Manopt_opts, 'solver'))
    % Use the trust-region solver by default
    Manopt_opts.solver = @trustregions;
end
solver_name = func2str(Manopt_opts.solver);
if (~strcmp(solver_name, 'trustregions') && ~strcmp(solver_name, 'conjugategradient') && ~strcmp(solver_name, 'steepestdescent'))
    error(sprintf('Unrecognized Manopt solver: %s', solver_name));
end
fprintf('\nSolving Riemannian optimization problems using Manopt''s "%s" solver\n\n', solver_name);

% Set cost function handles
manopt_data.cost = @(Y) evaluate_objective(Y', problem_data, SE_Sync_opts.Cholesky);
manopt_data.egrad = @(Y) Euclidean_gradient(Y', problem_data, SE_Sync_opts.Cholesky)';
manopt_data.ehess = @(Y, Ydot) Euclidean_Hessian_vector_product(Y', Ydot', problem_data, SE_Sync_opts.Cholesky)';

% Set preconditioning function, if desired
if(exist('precon', 'var'))
    manopt_data.precon = @(x,u) manopt_data.M.proj(x, precon(u));
end


% Set additional stopping criterion for Manopt: stop if the relative
% decrease in function value between successive iterates drops below the
% threshold specified in SE_Sync_opts.relative_func_decrease_tol
if(strcmp(solver_name, 'trustregions'))
    Manopt_opts.stopfun = @(manopt_problem, x, info, last) relative_func_decrease_stopfun(manopt_problem, x, info, last, Manopt_opts.rel_func_tol);
end

% Log the sequence of iterates visited by the Riemannian Staircase
Manopt_opts.statsfun = @log_iterates;


% Counter to keep track of how many iterations of the Riemannian Staircase
% have been performed
iter = 0;

%%  RIEMANNIAN STAIRCASE
for r = SE_Sync_opts.r0 : SE_Sync_opts.rmax
    iter = iter + 1;  % Increment iteration number
    
    % Starting at Y0, use Manopt's truncated-Newton trust-region method to
    % descend to a first-order critical point.
    
    fprintf('\nRIEMANNIAN STAIRCASE (level r = %d):\n', r);
    
    [YoptT, Fval, manopt_info, Manopt_opts] = manoptsolve(manopt_data, Y0', Manopt_opts);
    Yopt = YoptT';
    SDPLRval = Fval(end);
    
    % Store the optimal value and the elapsed computation time
    SDPLRvals(iter) = SDPLRval;
    optimization_times(iter) = manopt_info(end).time;
    
    % Store gradient norm and state traces
    gradnorms = [gradnorms, manopt_info.gradnorm];
    Yvals = [Yvals, {manopt_info.Yvals}];
    
    % Augment Yopt by padding with an additional row of zeros; this
    % preserves Yopt's first-order criticality while ensuring that it is
    % rank-deficient
    
    Yplus = vertcat(Yopt, zeros(1, problem_data.d * problem_data.n));
    
    
    fprintf('\nChecking second-order optimality...\n');
    % At this point, Yplus is a rank-deficient critial point, so check
    % 2nd-order optimality conditions
    
    % Compute Lagrange multiplier matrix Lambda corresponding to Yplus
    Lambda = compute_Lambda(Yopt, problem_data, SE_Sync_opts.Cholesky);
    
    % Compute minimum eigenvalue/eigenvector pair for Q - Lambda
    tic();
    [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, SE_Sync_opts.min_eig_num_tol, SE_Sync_opts.eig_comp_max_iters, SE_Sync_opts.Cholesky);
    min_eig_comp_time = toc();
    
    % Store the minimum eigenvalue and elapsed computation times
    min_eig_vals(iter) = lambda_min;
    min_eig_times(iter) = min_eig_comp_time;
    
    if( lambda_min > -SE_Sync_opts.min_eig_num_tol)
        % Yopt is a second-order critical point
        fprintf('Found second-order critical point! (minimum eigenvalue = %g, elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
        break;
    else
        fprintf('Saddle point detected (minimum eigenvalue = %g,  elapsed computation time %g seconds)\n', lambda_min, min_eig_comp_time);
        % lambda_min is a negative eigenvalue of Q - Lambda, so the KKT
        % conditions for the semidefinite relaxation are not satisfied;
        % this implies that Yplus is a saddle point of the rank-restricted
        % semidefinite optimization.  Fortunately, the eigenvector v
        % corresponding to lambda_min can be used to provide a descent
        % direction from this saddle point, as described in Theorem 3.9 of
        % the paper "A Riemannian Low-Rank Method for Optimization over
        % Semidefinite Matrices with Block-Diagonal Constraints".
        
        % Define the vector Ydot := e_{r+1} * v'; this is tangent to the
        % manifold St(d, r+1)^n at Yplus and provides a direction of
        % negative curvature
        disp('Computing escape direction...');
        Ydot = vertcat(zeros(r, problem_data.d * problem_data.n), v');
        
        % Augment the dimensionality of the Stiefel manifolds in
        % preparation for the next iteration
        
        manopt_data.M = stiefelstackedfactory(problem_data.n, problem_data.d, r+1);
        % Update preconditioning function, if it's used
        if(exist('precon', 'var'))
            manopt_data.precon = @(x,u) manopt_data.M.proj(x, precon(u));
        end
        
        % Perform line search along the escape direction Ydot to escape the
        % saddle point and obtain the initial iterate for the next level in
        % the Staircase
        
        % Compute a scaling factor alpha such that the scaled step
        % alpha*Ydot' should produce a trial point Ytest whose gradient has
        % a norm 100 times greater than the gradient tolerance stopping
        % criterion currently being used in the RTR optimization routine
        alpha = Manopt_opts.tolgradnorm / (norm(v) * abs(lambda_min));
        
        disp('Line searching along escape direction to escape saddle point...');
        tic();
        [stepsize, Y0T] = linesearch_decrease(manopt_data, Yplus', alpha * Ydot', SDPLRval);
        line_search_time = toc();
        Y0 = Y0T';
        fprintf('Line search completed (elapsed computation time %g seconds)\n', line_search_time);
    end
end

fprintf('\n\n===== END RIEMANNIAN STAIRCASE =====\n\n');

%% POST-PROCESSING

% Return optimal value of the SDP (in the case that a rank-deficient,
% second-order critical point is obtained, this is equal to the optimum
% value obtained from the Riemannian optimization

SDPval = SDPLRval;

disp('Rounding solution...');
% Round the solution
tic();
Rhat = round_solution(Yopt, problem_data);
solution_rounding_time = toc();
fprintf('Elapsed computation time: %g seconds\n\n', solution_rounding_time);

disp('Recovering translational estimates...');
% Recover the optimal translational estimates
tic();
that = recover_translations(Rhat, problem_data);
translation_recovery_time = toc();
fprintf('Elapsed computation time: %g seconds\n\n', translation_recovery_time);

xhat.R = Rhat;
xhat.t = that;

Fxhat = evaluate_objective(Rhat, problem_data, SE_Sync_opts.Cholesky);

fprintf('Value of SDP solution F(Y): %g\n', SDPval);
fprintf('Norm of Riemannian gradient grad F(Y): %g\n', manopt_info(end).gradnorm);
fprintf('Value of rounded pose estimate xhat: %g\n', Fxhat);
fprintf('Suboptimality bound of recovered pose estimate: %g\n\n', Fxhat - SDPval);
total_computation_time = toc(timerVal);

fprintf('Total elapsed computation time: %g seconds\n\n', total_computation_time);

% Output info
SE_Sync_info.mat_construct_times = auxiliary_matrix_construction_time;
SE_Sync_info.init_time = init_time;
SE_Sync_info.SDPLRvals = SDPLRvals(1:iter);
SE_Sync_info.optimization_times = optimization_times(1:iter);
SE_Sync_info.min_eig_vals = min_eig_vals(1:iter);
SE_Sync_info.min_eig_times = min_eig_times(1:iter);
SE_Sync_info.manopt_info = manopt_info;
SE_Sync_info.total_computation_time = total_computation_time;
SE_Sync_info.Yvals = Yvals;
SE_Sync_info.gradnorms = gradnorms;

fprintf('\n===== END SE-SYNC =====\n');


end