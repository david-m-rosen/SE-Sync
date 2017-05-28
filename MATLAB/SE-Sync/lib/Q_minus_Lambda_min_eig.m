function [lambda_min, v_min, flag] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, max_iters, use_Cholesky)
%function [lambda_min, v_min, flag] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, max_iters, use_Cholesky)
%
% Given the Lagrange multiplier Lambda corresponding to a critical point
% Yopt of the low-rank Riemannian optimization problem, this function
% computes and returns the minimum (algebraically smallest) eigenvalue
% of the matrix Q - Lambda, together with a corresponding eigenvector v.
% Here 'tol'  refers to the relative tolerance of the minimum eigenvalue computation
% using MATLAB's 'eigs' function.

% Copyright (C) 2016 by David M. Rosen

if nargin < 4
    tol = eps;  % default value
end

if(nargin >= 5)
    eigs_opts.maxit = max_iters;
end

if nargin < 6
    use_Cholesky = true;
end

% First, estimate the largest-magnitude eigenvalue of Q - Lambda
eigs_opts.issym = true;
eigs_opts.isreal = true;
eigs_opts.tol = tol;


% This function returns the product (Q - Lambda)*x
QminusLambda = @(x) Q_minus_Lambda_product(x, Lambda, problem_data, use_Cholesky);

[v_lm, lambda_lm, flag] = eigs(QminusLambda, problem_data.d*problem_data.n, 1, 'LM', eigs_opts);

if(flag ~= 0)
    warning('Minimum eigenvalue computation did not converge within the desired tolerance!\n');
end

if(lambda_lm < 0)
    % The largest-magnitude eigenvalue is negative, and therefore also the
    % minimum eigenvalue, so just return this solution
    
    lambda_min = lambda_lm;
    v_min = v_lm;
else
    % The largest-magnitude eigenvalue is positive, and therefore the
    % maximum eigenvalue.  Therefore, after shifting the spectrum of Q -
    % Lambda by -2*lambda_lm (by forming Q - Lambda - 2*lambda_max*I), the
    % shifted spectrum will line in the interval [lambda_min(A) - 2*
    % lambda_max(A), -lambda_max*A]; in particular, the largest-magnitude eigenvalue of
    % Q - Lambda - 2*lambda_max*I is lambda_min - 2*lambda_max, with
    % corresponding eigenvector v_min; furthermore, the condition number
    % sigma of Q - Lambda - 2*lambda_max is then upper-bounded by 2 :-).
    
    % Function to compute and return (Q - Lambda - 2*lambda_max*I) * x
    QminusLambda_shifted = @(x) QminusLambda(x) - 2*lambda_lm*x;
    
    if nargin >= 3
        % In the case that exactness holds, the minimum eigenvector
        % will be 0, with corresponding eigenvectors the columns of Yopt', so
        % we would like to use this as an initial guess.
        eigs_opts.v0 = Yopt(1, :)';
    end
    
    [v_min, shifted_lambda_min, flag] = eigs(QminusLambda_shifted, problem_data.d*problem_data.n, 1, 'LM', eigs_opts);
    lambda_min = shifted_lambda_min + 2*lambda_lm;
    if(flag ~= 0)
        warning('Minimum eigenvalue computation did not converge within the desired tolerance!\n');
    end
end
end

