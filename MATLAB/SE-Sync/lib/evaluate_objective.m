function [trQYtY, YQ] = evaluate_objective(Y, problem_data, use_Cholesky)
% function [trQYtY, YQ] = evaluate_objective(Y, problem_data, use_Cholesky)
%
% This function computes and returns the value of the objective function
% tr(Q Y^T Y).  Optionally, it returns the product YQ as the second
% argument

% Copyright (C) 2016 by David M. Rosen

if nargin < 3
    use_Cholesky = true;
end

Yt = Y';
YQ = Qproduct(Yt, problem_data, use_Cholesky)';

trQYtY = trace(YQ * Yt);


end

