function [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, problem_data, tol, use_Cholesky)
%
%function [lambda_min, v] = Q_minus_Lambda_min_eig(Lambda, problem_data, tol, use_Cholesky)
%
% Given the Lagrange multiplier Lambda corresponding to a critical point
% Yopt of the low-rank Riemannian optimization problem, this function
% computes and returns the minimum eigenvalue lambda_min of the matrix 
% Q - Lambda, together with an eigenvector v corresponding to this
% eigenvalue.  Here 'tol' refers to the relative tolerance of the minimum
% eigenvalue computation using MATLAB's 'eigs' function

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    tol = 1e-5;  % default value
end

if nargin < 4
    use_Cholesky = true;
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
eigs_opts.tol = tol;  %This should be a reasonably sharp estimate
[v, shifted_lambda_min] = eigs(QminusLambda_shifted, problem_data.d*problem_data.n, 1, 'SA', eigs_opts);
lambda_min = shifted_lambda_min - lambda_max;

end

