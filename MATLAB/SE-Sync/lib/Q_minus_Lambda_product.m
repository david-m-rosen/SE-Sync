function Q_minus_Lambda_X = Q_minus_Lambda_product(X, Lambda, problem_data, use_Cholesky)
%function Q_minus_Lambda_x = Q_minus_Lambda_product(X, Lambda, problem_data, use_Cholesky)
%
% This function computes and returns the product (Q - Lambda)*X.

% Copyright (C) 2016 by David M. Rosen

if nargin < 4
    use_Cholesky = true;
end

Lambda_X = zeros(problem_data.d*problem_data.n, size(X, 2));

for i = 1:problem_data.n
    rmin = problem_data.d*(i-1) + 1;
    rmax = problem_data.d*(i-1) + problem_data.d;
    Lambda_X(rmin:rmax, :) = Lambda(:, rmin:rmax) * X(rmin:rmax, :);
end

Q_minus_Lambda_X = Qproduct(X, problem_data, use_Cholesky) - Lambda_X;

end