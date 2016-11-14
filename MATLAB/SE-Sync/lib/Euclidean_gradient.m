function egrad = Euclidean_gradient(Y, problem_data, use_Cholesky)
%function egrad = Euclidean_gradient(Y, problem_data, use_Cholesky)
%
% This function computes and returns the value of the Euclidean gradient of
% the objective function: nabla F(Y) = 2YQ.

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    use_Cholesky = true;
end

egrad = 2*Qproduct(Y', problem_data, use_Cholesky)';



end

