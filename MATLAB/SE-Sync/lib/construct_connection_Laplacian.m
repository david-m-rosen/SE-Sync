function Lrho = construct_connection_Laplacian(measurements)
%function Lrho = construct_connection_Laplacian(measurements)
%
% This function computes and returns the connection Laplacian for the
% rotational observations (see eq. (15) of the paper).

D = length(measurements.t{1}); % D = dimension of SE(d)
N = max(max(measurements.edges));  % N = number of nodes in the pose graph
M = size(measurements.edges,1); % M = number of edges in the pose graph

% The number of nonzero elements in the connection Laplacian; there are 
% 2*D^2*M nonzero elements corresponding to the off-diagonal blocks 
% containing the raw measurements (one dxd rotation matrix above the main 
% diagonal, and another below, for each of the M measurements, plus D*N 
% nonzero elements along the main diagonal for the scaled identity
% matrices.  Note that in the case of multiple observations between two
% nodes x_i and x_j the corresponding (i,j) block in the connection
% Laplacian should contain the *sum* of the weighted observations
% kappa_{ij} R_{ij} between these nodes.  We achieve this result below by
% exploiting the fact that MATLAB's sparse() command forms the resulting
% sparse matrix by *summing* values with the same (i,j) indices.

% Copyright (C) 2016 by David M. Rosen

D2 = D^2;
off_diag_inc = 2*D2;

NNZ = off_diag_inc * M + D*N;

%Allocate storage for the row, column, and value vectors

rows = zeros(1, NNZ);
cols = zeros(1, NNZ);
vals = zeros(1, NNZ);

degs = zeros(1, N);  %Vector to store the degrees of the nodes

%Iterate over the measurements in the pose graph
for k = 1:M
   
    %EXTRACT MEASUREMENT DATA
    i = measurements.edges(k, 1);  %The node that this edge leaves
    j = measurements.edges(k, 2);  %The node that this edge enters
    
    Rij = measurements.R{k};  %The rotation matrix for this observation
    
    kappa = measurements.kappa{k};  %The precision for this rotational observation
    
    
    %PROCESS MEASUREMENT DATA
    
    %Increment the degrees of the ith and jth nodes by kappa
    degs(i) = degs(i) + kappa;
    degs(j) = degs(j) + kappa;
    
    %Set the (i,j)th DxD block of Lrho = -kappa*Rij
    [r, c, Rvect] = rcvize_matrix(Rij, i, j);
    rows(off_diag_inc*(k - 1) + 1 : off_diag_inc*(k-1) + D2) = r;
    cols(off_diag_inc*(k - 1) + 1 : off_diag_inc*(k-1) + D2) = c;
    vals(off_diag_inc*(k - 1) + 1 : off_diag_inc*(k-1) + D2) = -kappa*Rvect;
    
    %Set the (j,i)th DxD block of Lrho = -kappa*R_ij^T
    [r, c, Rvect] = rcvize_matrix(Rij', j, i);
    rows(off_diag_inc*(k - 1) + D2 + 1 : off_diag_inc*(k-1) + off_diag_inc) = r;
    cols(off_diag_inc*(k - 1) + D2 + 1 : off_diag_inc*(k-1) + off_diag_inc) = c;
    vals(off_diag_inc*(k - 1) + D2 + 1 : off_diag_inc*(k-1) + off_diag_inc) = -kappa*Rvect;
end

%Now set the diagonal elements to be D-fold copies of the degrees of each
%node

rows(off_diag_inc*M + 1 : NNZ) = [1:D*N];
cols(off_diag_inc*M + 1 : NNZ) = [1:D*N];
vals(off_diag_inc*M + 1 : NNZ) = kron(degs, ones(1,D));

Lrho = sparse(rows, cols, vals, D*N, D*N);


end

