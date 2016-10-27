function egrad = Euclidean_gradient(Y, problem_data)
%function egrad = Euclidean_gradient(Y, problem_data)
%
% This function computes and returns the value of the Euclidean gradient of
% the objective function: nabla F(Y) = 2YQ.

egrad = 2*Qproduct(Y', problem_data)';



end

