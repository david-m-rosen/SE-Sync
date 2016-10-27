function grad = getGradient(problem, x, storedb, key)
% Computes the gradient of the cost function at x.
%
% function grad = getGradient(problem, x)
% function grad = getGradient(problem, x, storedb)
% function grad = getGradient(problem, x, storedb, key)
%
% Returns the gradient at x of the cost function described in the problem
% structure.
%
% storedb is a StoreDB object, key is the StoreDB key to point x.
%
% See also: getDirectionalDerivative canGetGradient

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 
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

    
    if isfield(problem, 'grad')
    %% Compute the gradient using grad.
	
        % Check whether this function wants to deal with storedb or not.
        switch nargin(problem.cost)
            case 1
                grad = problem.grad(x);
            case 2
                % Obtain, pass along, and save the store for x.
                store = storedb.getWithShared(key);
                [grad, store] = problem.grad(x, store);
                storedb.setWithShared(store, key);
            case 3
                % Pass along the whole storedb (by reference), with key.
                grad = problem.grad(x, storedb, key);
            otherwise
                up = MException('manopt:getGradient:badgrad', ...
                    'grad should accept 1, 2 or 3 inputs.');
                throw(up);
        end
    
    elseif isfield(problem, 'costgrad')
    %% Compute the gradient using costgrad.
		
        % Check whether this function wants to deal with storedb or not.
        switch nargin(problem.costgrad)
            case 1
                [unused, grad] = problem.costgrad(x); %#ok
            case 2
                % Obtain, pass along, and save the store for x.
                store = storedb.getWithShared(key);
                [unused, grad, store] = problem.costgrad(x, store); %#ok
                storedb.setWithShared(store, key);
            case 3
                % Pass along the whole storedb (by reference), with key.
                [unused, grad] = problem.costgrad(x, storedb, key); %#ok
            otherwise
                up = MException('manopt:getGradient:badcostgrad', ...
                    'costgrad should accept 1, 2 or 3 inputs.');
                throw(up);
        end
    
    elseif canGetEuclideanGradient(problem)
    %% Compute the gradient using the Euclidean gradient.
        
        egrad = getEuclideanGradient(problem, x, storedb, key);
        grad = problem.M.egrad2rgrad(x, egrad);

    else
    %% Abandon computing the gradient.
    
        up = MException('manopt:getGradient:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute the gradient of the cost.']);
        throw(up);
        
    end
    
end
