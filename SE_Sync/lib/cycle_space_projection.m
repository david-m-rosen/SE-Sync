function Mproj = cycle_space_projection(M, Ared, L)
%function Mproj = cycle_space_projection(M, Ared, L)
%
%This function computes the product Pi*M, where Pi is the orthogonal
%projection matrix onto the cycle space of a graph G, given:
%
% -Ared:  A reduced oriented incidence matrix for the graph, and 
% -L: A sparse lower-triangular Cholesky factor for the symmetric 
% positive-definite product (Ared*Ared')

P1 = Ared*M;
P2 = L \ P1;
P3 = L' \ P2;
P4 = Ared' * P3;

Mproj = M - P4;

end

