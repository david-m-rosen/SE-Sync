function Hvec = Euclidean_Hessian_vector_product(Y, Ydot, problem_data, use_Cholesky)
%function Hvec = Euclidean_Hessian_vector_product(Y, Ydot, problem_data, use_Cholesky)
%
% This function computes and returns the value of the Euclidean Hessian at
% the point Y evaluated along the tangent direction Ydot.

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    use_Cholesky = true;
end

Hvec = 2*Qproduct(Ydot', problem_data, use_Cholesky)';



end

