function QX = Qproduct(X, problem_data)
%function QX = Qproduct(X, problem_data)
%
%This function computes and returns the matrix product P = Q*M, where 
%
%Q = ConLap + T'* Pi * diag(Omega) * Pi * T
%
%is the quadratic form that defines the objective function for the
%pose-graph optimization problem.  The arguments to this function are:
%
% -ConLap:  The (sparse) connection Laplacian for the subgraph of rotational
% observations
% -Ared:  The (sparse) reduced reduced adjacency matrix for the underlying
% (unweighted) pose graph
% -L:  The (sparse) lower-triangular Cholesky factor for the reduced
% (unweighted) graph Laplacian LapRed = (Ared * Ared')
% -T:  The (sparse) matrix constructed from the raw translational
% observations
% -Omega: The (sparse) matrix whose main diagonal contains the precisions
% of the translational measurements


%TRANSLATIONAL TERM

%We begin by computing the translational term:
%
% Qtau = T' * Pi * Omega * Pi * T * M
%
% In order to avoid having to store and/or manipulate prohibitively large
% dense matrices, we compute this product by working associatively from
% right to left

P1 = problem_data.T*X;  
P2 = cycle_space_projection(P1, problem_data.Ared, problem_data.L);
P3 = problem_data.Omega * P2;
P4 = cycle_space_projection(P3, problem_data.Ared, problem_data.L);
Qtau = problem_data.T' * P4;

%Now multiply the rows of

QX = problem_data.ConLap*X + Qtau;

end

