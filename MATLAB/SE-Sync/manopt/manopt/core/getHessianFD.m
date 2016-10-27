function hessfd = getHessianFD(problem, x, d, storedb, key)
% Computes an approx. of the Hessian w/ finite differences of the gradient.
%
% function hessfd = getHessianFD(problem, x, d)
% function hessfd = getHessianFD(problem, x, d, storedb)
% function hessfd = getHessianFD(problem, x, d, storedb, key)
%
% Returns a finite difference approximation of the Hessian at x along d of
% the cost function described in the problem structure. The finite
% difference is based on computations of the gradient.
%
% storedb is a StoreDB object, key is the StoreDB key to point x.
%
% If the gradient cannot be computed, an exception is thrown.
%
% See also: approxhessianFD

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
%
%   Feb. 19, 2015 (NB):
%       It is sufficient to ensure positive radial linearity to guarantee
%       (together with other assumptions) that this approximation of the
%       Hessian will confer global convergence to the trust-regions method.
%       Formerly, in-code comments referred to the necessity of having
%       complete radial linearity, and that this was harder to achieve.
%       This appears not to be necessary after all, which simplifies the
%       code.
%
%   April 3, 2015 (NB):
%       Works with the new StoreDB class system.

    % Allow omission of the key, and even of storedb.
    if ~exist('key', 'var')
        if ~exist('storedb', 'var')
            storedb = StoreDB();
        end
        key = storedb.getNewKey();
    end

    
    if ~canGetGradient(problem)
        up = MException('manopt:getHessianFD:nogradient', ...
            'getHessianFD requires the gradient to be computable.');
        throw(up);
    end
	
	% Step size
    norm_d = problem.M.norm(x, d);
    
    % First, check whether the step d is not too small
    if norm_d < eps
        hessfd = problem.M.zerovec(x);
        return;
    end
    
    % Parameter: how far do we look?
    % (Use approxhessianFD explicitly to gain access to this parameter.)
    epsilon = 1e-4;
        
    c = epsilon/norm_d;
    
    % Compute the gradient at the current point.
    grad = getGradient(problem, x, storedb, key);
    
    % Compute a point a little further along d and the gradient there.
    % Since this is a new point, we need a new key for it, for the storedb.
    x1 = problem.M.retr(x, d, c);
    key1 = storedb.getNewKey();
    grad1 = getGradient(problem, x1, storedb, key1);
    
    % Transport grad1 back from x1 to x.
    grad1 = problem.M.transp(x1, x, grad1);
    
    % Return the finite difference of them.
    hessfd = problem.M.lincomb(x, 1/c, grad1, -1/c, grad);
    
end
