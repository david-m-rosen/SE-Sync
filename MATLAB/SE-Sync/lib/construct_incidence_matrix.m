function A = construct_incidence_matrix(measurements)
%function A = construct_incidence_matrix(measurements)
%
% This function computes and returns the oriented incidence matrix of the
% underlying directed graph of measurements:
%
% A_{ie} =  -1, if edge e *leaves* node i,
%           +1, if edge e *enters* node i,
%           0, otherwise.
%
% (see eq. (7) in the paper).

% Copyright (C) 2016 by David M. Rosen

N = max(max(measurements.edges));  % Number of nodes in the pose graph
M = size(measurements.edges, 1);  % Number of edges in the pose graph

out_nodes = measurements.edges(:, 1)';  %out_nodes(e) = i if edge e leaves node i
in_nodes = measurements.edges(:, 2)';  %out_nodes(e) = j if edge e enters node i

node_indices = [out_nodes, in_nodes];
edge_indices = [ [1:M], [1:M] ];
vals = [-ones(1, M), ones(1, M)];

A = sparse(node_indices, edge_indices, vals, N, M);


end

