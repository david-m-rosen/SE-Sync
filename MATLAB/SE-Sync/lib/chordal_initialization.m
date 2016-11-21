function [Rchordal, tchordal] = chordal_initialization(measurements)
% function [Rchordal, tchordal] = chordal_initialization(measurements)
%
% This helper function accepts a MATLAB struct containing the raw 
% measurements specifying a special Euclidean synchronization problem, and 
% returns the corresponding chordal initialization for this problem.
%
% INPUT: A MATLAB struct 'measurements' containing the following fields
% (see eq. (11) in the long-form version of the paper for details):
% edges:  An (mx2)-dimension encoding the edges in the measurement network;
%     edges(k, :) = [i,j] means that the kth measurement is of the
%     relative transform from pose i to pose j.  NB:  This indexing scheme 
%     requires that the states x_i are numbered sequentially as 
%     x_1, ... x_n.
% R:  An m-dimensional cell array whose kth element is the rotational part
%     of the kth measurement
% t:  An m-dimensional cell array whose kth element is the translational
%     part of the kth measurement
% kappa:  An m-dimensional cell array whose kth element gives the precision
%     of the rotational part of the kth measurement. 
% tau:  An m-dimensional cell array whose kth element gives the precision
%     of the translational part of the kth measurement.

% OUTPUT: 
%  -Rchordal:  A d x dn block matrix Rchordal containing the estimating
%     orientations.
%  -tchordal [optional]:  A d x n block matrix containing the estimated
%     positions.

% Copyright (C) 2016 by David M. Rosen

d = length(measurements.t{1});
n = max(max(measurements.edges));
m = size(measurements.edges, 1);

if nargout > 1
    [B3, B2, B1] = construct_B_matrices(measurements);
else
    B3 = construct_B_matrices(measurements);
end


% First, estimate the rotations using *only* the rotational observations
Id = eye(d);
Id_vec = Id(:);

% Compute the constant vector cR induced by fixing the first orientation
% estimate to be the identity Id.
cR = B3(:, 1:d^2) * Id_vec;

% Compute an estimate of the remaining rotations by solving the resulting
% least squares problem *without* enforcing the constraint the the
% estimates lie in SO(d)
r2vec = - B3(:, d^2 + 1 : end) \ cR;
rvec = vertcat(Id_vec, r2vec);
R_LS = reshape(rvec, d, d*n);

% Now reproject these estimates onto SO(d)
Rchordal = zeros(d, d*n);
for i = 1:n
    Rchordal(:, d*(i-1) + 1 : d*i) = project_to_SOd(R_LS(:, d*(i-1) + 1 : d*i));
end

if nargout > 1
    
    % Solve for the translations in terms of the rotations
    
    % Constant vector induced by fixing R
    cT = B2 * Rchordal(:);
    
    % Solve for t_2, ... t_n assuming that t_1 = 0.
    
    t2 = - B1(:, d + 1 : end) \ cT;
    tvec = vertcat(zeros(d, 1), t2);
    tchordal = reshape(tvec, d, n);
end

end

