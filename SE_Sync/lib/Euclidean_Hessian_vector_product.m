function Hvec = Euclidean_Hessian_vector_product(Y, Ydot, problem_data)
%function Hvec = Euclidean_Hessian_vector_product(Y, Ydot, problem_data)
%
% This function computes and returns the value of the Euclidean Hessian at
% the point Y evaluated along the tangent direction Ydot.

Hvec = 2*Qproduct(Ydot', problem_data)';



end

