function [lambdas, V] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, use_Cholesky, num_eigs)
%function [lambdas, V] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, use_Cholesky, num_eigs)
%
% Given the Lagrange multiplier Lambda corresponding to a critical point
% Yopt of the low-rank Riemannian optimization problem, this function
% computes and returns the num_eigs smallest eigenvalues of the matrix 
% Q - Lambda, together with their corresponding eigenvectors V. Here 'tol' 
% refers to the relative tolerance of the minimum eigenvalue computation 
% using MATLAB's 'eigs' function.

% Copyright (C) 2016 by David M. Rosen

if nargin < 4
    tol = 1e-5;  % default value
end

if nargin < 5
    use_Cholesky = true;
end

if nargin < 6
    num_eigs = 1;
end

% First, estimate the maximum eigenvalue of Q - Lambda (this should be its
% norm in the typical case)
eigs_opts.issym = true;
eigs_opts.isreal = true;

% This function returns the product (Q - Lambda)*x
QminusLambda = @(x) Q_minus_Lambda_product(x, Lambda, problem_data, use_Cholesky);

eigs_opts.tol = 100*tol;  %This estimate needn't be particularly sharp...
lambda_max = eigs(QminusLambda, problem_data.d*problem_data.n, 1, 'LA', eigs_opts);

% We shift the spectrum of Q - Lambda by adding lambda_max_est*I; this
% improves the condition number (to ~2 in the typical case)
% *without* perturbing any of the eigenspaces in Q - Lambda; this has the
% effect of producing MUCH faster convergence when running the Lanczos
% algorithm
%
QminusLambda_shifted = @(x) QminusLambda(x) + lambda_max*x;

% Now compute the minimum eigenvalue of Q - Lambda + lambda_max_est * I
eigs_opts.tol = tol;  %T his should be a reasonably sharp estimate

if nargin >= 3
    % In the (typical) case that exactness holds, the minimum eigenvector
    % will be 0, with corresponding eigenvectors the columns of Yopt', so
    % we would like to use a guess close to this.  However, we would not
    % like to use /exactly/ this guess, because we know that (Q - Lambda)Y'
    % = 0 implies that Y' lies in the null space of Q - Lambda, and
    % therefore an iterative algorithm will get "stuck" if we start
    % /exactly/ there.  Therefore, we will "fuzz" this initial guess by
    % adding a randomly-sampled perturbation that is small in norm relative
    % to the first column of Yopt; this enables us to excite modes other
    % than Yopt itself (thereby escaping from this subspace in the 'bad
    % case' that Yopt is not the minimum eigenvalue subspace), while still
    % remaining close enough that we can converge to this answer quickly in
    % the 'good' case
    
    relative_perturbation = .03;
    eigs_opts.v0 = Yopt(1, :)' + (relative_perturbation / sqrt(problem_data.d))*randn(problem_data.n * problem_data.d, 1);
end

[V, shifted_lambda_min] = eigs(QminusLambda_shifted, problem_data.d*problem_data.n, num_eigs, 'SA', eigs_opts);
lambdas = shifted_lambda_min - lambda_max*eye(num_eigs);

end

