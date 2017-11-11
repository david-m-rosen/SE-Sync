function [lambda_min, v_min, flag] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, max_iters, use_Cholesky)
%function [lambda_min, v_min, flag] = Q_minus_Lambda_min_eig(Lambda, problem_data, Yopt, tol, max_iters, use_Cholesky)
%
% Given the Lagrange multiplier Lambda corresponding to a critical point
% Yopt of the low-rank Riemannian optimization problem, this function
% computes and returns the minimum (algebraically smallest) eigenvalue
% of the matrix Q - Lambda, together with a corresponding eigenvector v.
% Here 'tol' refers to the (absolute!) tolerance of the minimum eigenvalue computation

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
eigs_opts.tol = 1e-3;  % Find largest eigenvalue to an accuracy of .1%


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
        % If Ystar is a critical point of F, then Ystar^T is also in the 
        % null space of Q - Lambda(Ystar) (cf. Lemma 6 of the tech report), 
        % and therefore its rows are eigenvectors corresponding to the 
        % eigenvalue 0.  In the case that the relaxation is exact, this is 
        % the *minimum* eigenvalue, and therefore the rows of Ystar are 
        % exactly the eigenvectors that we're looking for.  On the other 
        % hand, if the relaxation is *not* exact, then Q - Lambda(Ystar)
        % has at least one strictly negative eigenvalue, and the rows of 
        % Ystar are *unstable fixed points* for the Lanczos iterations.  
        % Thus, we will take a slightly "fuzzed" version of the first row 
        % of Ystar as an initialization for the Lanczos iterations; this 
        % allows for rapid convergence in the case that the relaxation is 
        % exact (since are starting close to a solution), while 
        % simultaneously allowing the iterations to escape from this fixed 
        % point in the case that the relaxation is not exact.
        
        v = Yopt(1, :)';
        
        % Determine the standard deviation necessary such that a vector of
        % n*d iid elements sampled from N(0, sigma^2) will perturb v by ~3%
        sigma = .03 * norm(v) / (problem_data.n * problem_data.d); 
        eigs_opts.v0 = v + sigma*randn(problem_data.n * problem_data.d, 1);
    end
    
    % Find minimum eigenvalue to within an absolute tolerance of 'tol'
    eigs_opts.tol = tol / lambda_lm;
    
    [v_min, shifted_lambda_min, flag] = eigs(QminusLambda_shifted, problem_data.d*problem_data.n, 1, 'LM', eigs_opts);
    lambda_min = shifted_lambda_min + 2*lambda_lm;
    if(flag ~= 0)
        warning('Minimum eigenvalue computation did not converge within the desired tolerance!\n');
    end
end
end

