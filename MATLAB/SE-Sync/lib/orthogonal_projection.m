function PiX = orthogonal_projection(X, problem_data, Cholesky )
% function PiX = orthogonal_projection(X, problem_data, Cholesky )
% 
% This function computes and returns the orthogonal projection of X onto
% ker(A * Omega^(1/2), using either the Cholesky factor L for the reduced
% Laplacian L(W^tau) or by applying an orthogonal (QR) decomposition for
% Omega^(1/2) * Ared'

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    Cholesky = true;
end


if Cholesky
    % Compute the projection using a sequence of matrix multiplications and
    % upper-triangular solves
    
    P1 = problem_data.sqrt_Omega_AredT' * X;
    P2 = problem_data.L \ P1;
    P3 = problem_data.L' \ P2;
    P4 = problem_data.sqrt_Omega_AredT * P3;
    
    PiX = X - P4;
    
else
    % Use QR decomposition
    
    vstar = problem_data.sqrt_Omega_AredT \ X;
    PiX = X - problem_data.sqrt_Omega_AredT * vstar;
end


end

