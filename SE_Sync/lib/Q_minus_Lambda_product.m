function Q_minus_Lambda_x = Q_minus_Lambda_product(x, Lambda, problem_data )
%function Q_minus_Lambda_x = Q_minus_Lambda_product(x, Lambda, problem_data )
%
% This function computes and returns the matrix-vector product (Q - Lambda)*x.

Lambda_x = zeros(size(x));

for i = 1:problem_data.n
    cmin = problem_data.d*(i-1) + 1;
    cmax = problem_data.d*(i-1) + problem_data.d;
    Lambda_x(cmin:cmax) = Lambda(:, cmin:cmax) * x(cmin:cmax);
end

Q_minus_Lambda_x = Qproduct(x, problem_data) - Lambda_x;

end

