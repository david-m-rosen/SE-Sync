function [rows, cols, vals] = rcvize_matrix(M, i, j)
%function [rows, cols, vals] = rcvize_matrix(M, i, j)
%
% Given a DxD matrix M that we wish to insert into the (i,j)th DxD
% block of a sparse matrix, this function computes and returns the
% corresponding {row, col, val} vectors describing this block

% Copyright (C) 2016 by David M. Rosen

D = size(M, 1);

vals = reshape(M, [1, D^2]);  %Vectorize M by concatenating its columns

rows = repmat( [D*(i-1) + 1 : D*(i-1) + D], [1,D]);
cols = kron( [D*(j-1) + 1 : D*(j-1) + D], ones(1, D));


end

