function P = cycle_space_projection(M, Ared, L, Cholesky)
%function P = cycle_space_projection(M, Ared, L, Cholesky)
%
%This function computes the product Pi*M, where Pi is the orthogonal
%projection matrix onto the cycle space of a graph G.  
%
% Arguments:
%
% -Ared:  A reduced oriented incidence matrix for the graph G. 
%
% -L [optional]: A sparse lower-triangular Cholesky factor for the symmetric 
% positive-definite product (Ared*Ared')
% 
% -Cholesky [optional]:  A Boolean value indicating whether to compute the
% projection P = Pi * M via backsubstitution using the (cached) Cholesky 
% factor L, or via an orthogonal (QR) decomposition.  The former is often
% faster, but the latter is more numerically stable.
%
% Note that both the L and Cholesky arguments are optional.

% Copyright (C) 2016 by David M. Rosen

if( (nargin < 3) | ((nargin == 4 & ~Cholesky)))
    % Use QR decomposition to compute product X = (Ared * Ared')^{-1} * Ared * M
    %
    % Note that since Ared is full-rank, this corresponds to the solution
    % of the linear least-squares problem:
    %
    % || Ared' X - M ||_2
    %
    % via the normal equations
    
    X = Ared' \ M;
else
    % Compute the product
    %
    % X = (Ared * Ared')^{-1} * Ared * M
    %  = L^{-T} * L^{-1} * Ared * M
    %
    % via backsubstitution using the cached Cholesky factor L
    
    Y1 = Ared*M;
    Y2 = L \ Y1;
    X = L' \ Y2;
end

P = M - Ared' * X;

end

