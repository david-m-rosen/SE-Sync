function egrad = getEuclideanGradient(problem, x, storedb, key)
% Computes the Euclidean gradient of the cost function at x.
%
% function egrad = getEuclideanGradient(problem, x)
% function egrad = getEuclideanGradient(problem, x, storedb)
% function egrad = getEuclideanGradient(problem, x, storedb, key)
%
% Returns the Euclidean gradient at x of the cost function described in the
% problem structure.
%
% storedb is a StoreDB object, key is the StoreDB key to point x.
%
% Because computing the Hessian based on the Euclidean Hessian will require
% the Euclidean gradient every time, to avoid overly redundant
% computations, if the egrad function does not use the store caching
% capabilites, this implements an automatic caching functionality. Writing
% egrad to accept the optional store or storedb parameter will disable
% automatic caching, but allow user controlled caching.
%
% See also: getGradient canGetGradient canGetEuclideanGradient

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 9, 2013.
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

    
    if isfield(problem, 'egrad')
    %% Compute the Euclidean gradient using egrad.
	
        % Check whether this function wants to deal with storedb or not.
        switch nargin(problem.egrad)
            case 1
                % If it does not want to deal with the store structure,
                % then we do some caching of our own. There is a small
                % performance hit for this is some cases, but we expect
                % that this is most often the preferred choice.
                store = storedb.get(key);
                if ~isfield(store, 'egrad__')
                    store.egrad__ = problem.egrad(x);
                    storedb.set(store, key);
                end
                egrad = store.egrad__;
            case 2
                % Obtain, pass along, and save the store for x.
                % If the user deals with the store structure, then we don't
                % do any automatic caching: the user is in control.
                store = storedb.getWithShared(key);
                [egrad, store] = problem.egrad(x, store);
                storedb.setWithShared(store, key);
            case 3
                % Pass along the whole storedb (by reference), with key.
                % Same here: no automatic caching.
                egrad = problem.egrad(x, storedb, key);
            otherwise
                up = MException('manopt:getEuclideanGradient:badegrad', ...
                    'egrad should accept 1, 2 or 3 inputs.');
                throw(up);
        end

    else
    %% Abandon computing the Euclidean gradient
    
        up = MException('manopt:getEuclideanGradient:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute the Euclidean gradient of the cost.']);
        throw(up);
        
    end
    
end
