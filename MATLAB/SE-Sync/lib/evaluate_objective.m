function [trQYtY, YQ] = evaluate_objective(Y, problem_data)
% function [trQYtY, YQ] = evaluate_objective(Y, problem_data)
%
% This function computes and returns the value of the objective function
% tr(Q Y^T Y).  Optionally, it returns the product YQ as the second
% argument

% Copyright (C) 2016 by David M. Rosen

Yt = Y';
YQ = Qproduct(Yt, problem_data)';

trQYtY = trace(YQ * Yt);


end

