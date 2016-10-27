function costs = surfprofile(problem, x, d1, d2, t1, t2)
% Plot the cost function as a surface over a 2-dimensional subspace.
%
% function surfprofile(problem, x, d1, d2, t1, t2)
% function costs = surfprofile(problem, x, d1, d2, t1, t2)
%
% Evaluates the cost function at points
%
%   gamma(t1, t2) = exponential_x(t1*d1 + t2*d2)
% 
% where the exponential map at x is specified by problem.M.exp (retr is
% used instead if needed). d1 and d2 are two tangent vectors to problem.M
% at the point x. The values assigned to t1 and t2 are as specified in the
% two input vectors t1 and t2.
% 
% If the function is called with an output, the plot is not drawn and the
% values of the cost are returned in a matrix of size
% length(t1)*length(t2). To plot a surf, call surf(t1, t2, costs.') (notice
% the transpose).

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Sep. 1, 2014.
% Contributors: 
% Change log: 
%
%   April 3, 2015 (NB):
%       Works with the new StoreDB class system.

    % Verify that the problem description is sufficient.
    if ~canGetCost(problem)
        error('It seems no cost was provided.');  
    end
    
    if isfield(problem.M, 'exp')
        expo = problem.M.exp;
        str = 'Exp';
    else
        expo = problem.M.retr;
        str = 'Retr';
    end
    
    storedb = StoreDB();
    linesearch_fun = @(ta, tb) getCost(problem, ...
                         expo(x, problem.M.lincomb(x, ta, d1, tb, d2)), ...
                         storedb);
    
    costs = zeros(length(t1), length(t2));
    for i = 1 : length(t1)
        for j = 1 : length(t2)
            costs(i, j) = linesearch_fun(t1(i), t2(j));
        end
    end
    
    if nargout == 0
        surf(t1, t2, costs.');
        xlabel('t1');
        ylabel('t2');
        zlabel(['f(' str '_x(t1*d1+t2*d2))']);
    end
    
end
