function QX = Qproduct(X, problem_data, use_Cholesky)
% function QX = Qproduct(X, problem_data, use_Cholesky)
%
% This function computes and returns the matrix product P = Q*X, where 
%
% Q = L(G^rho) + Q^tau
%
% is the quadratic form that defines the objective function for the
% pose-graph optimization problem.

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    use_Cholesky = true;
end

% We begin by computing the translational term:
%
% Qtau = T' * Omega^(1/2) * Pi * Omega^(1/2) * T * X
%
% In order to avoid having to store and/or manipulate prohibitively large
% dense matrices, we compute this product by working associatively from
% right to left

QtauX = problem_data.sqrt_Omega_T' * orthogonal_projection(problem_data.sqrt_Omega_T * X, problem_data, use_Cholesky);


QX = problem_data.ConLap*X + QtauX;

end

