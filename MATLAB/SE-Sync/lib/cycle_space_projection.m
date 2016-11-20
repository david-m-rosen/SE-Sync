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

% Instead of computing P directly, we will instead compute its orthogonal
% complement, using the fact that ker(A)^\perp = image(A^T).  Therefore,
% the orthogonal complement of P is the vector realizing the minimum
% distance in:
%
% min_y || A^T y - x ||_2

Y = Ared' \ M;

V = Ared' * Y;

P = M - V;

end

