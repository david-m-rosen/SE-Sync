classdef StoreDB < handle_light
% The StoreDB class is a handle class to manage caching in Manopt.
%
% To create an object, call: storedb = StoreDB();
% Alternatively, call: storedb = StoreDB(storedepth); to instruct
% the database to keep at most storedepth store's in its history.
% (Note that clean up only happens when purge() is called).
%
% The storedb object is passed by reference: when it is passed to a
% function as an input, and that function modifies it, the original
% object is modified.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, April 3, 2015.
% Contributors: 
% Change log: 

% TODO : protect get/setWithShared calls: limit to one, and forbid access
%        to shared memory while it has not been returned.
%        Do think of the applyStatsFun case : calls a getWithShared, does
%        not need a setWithShared. I think for statsfun there should be a
%        method "forfeitWithShared".
    
    properties(Access = public)
       
        % This memory is meant to be shared at all times. Users can modify
        % this at will. It is the same for all points x.
        shared = struct();
        
        % This memory is used by the toolbox for, e.g., automatic caching
        % and book keeping. Users should not overwrite this. It is the
        % same for all points x.
        internal = struct();
        
        % When calling purge(), only a certain number of stores will be
        % kept in 'history'. This parameter fixes that number. The most
        % recently modified stores are kept. Set to inf to keep all stores.
        storedepth = inf;
        
    end
    
    properties(Access = private)
        
        % This structure holds separate memories for individual points.
        % Use get and set to interact with this. The field name 'shared' is
        % reserved, for use with get/setWithShared.
        history = struct();
        
        % This internal counter is used to obtain unique key's for points.
        counter = uint32(0);
        
        % This internal counter is used to time calls to 'set', and hence
        % keep track of which stores in 'history' were last updated.
        timer = uint32(0);
        
    end
    
    
    methods(Access = public)
        
        % Constructor
        function storedb = StoreDB(storedepth)
            if nargin >= 1
                storedb.storedepth = storedepth;
            end
        end
        
        % Return the store associated to a given key.
        % If the key is unknown, returns an empty structure.
        function store = get(storedb, key)
            if isfield(storedb.history, key)
                store = storedb.history.(key);
            else
                store = struct();
            end
        end
        
        % Same as get, but adds the shared memory in store.shared.
        function store = getWithShared(storedb, key)
            store = storedb.get(key);
            store.shared = storedb.shared;
        end
        
        % Save the given store at the given key. If no key is provided, a
        % new key is generated for this store (i.e., it is assumed this
        % store pertains to a new point). The key is returned in all cases.
        % A field 'lastset__' is added/updated in the store structure,
        % keeping track of the last time that store was modified.
        function key = set(storedb, store, key)
            if nargin < 3
                key = getNewKey(storedb);
            end
            store.lastset__ = storedb.timer;
            storedb.timer = storedb.timer + 1;
            storedb.history.(key) = store;
        end
        
        % Same as set, but extracts the shared memory and saves it.
        % The stored store will still have a 'shared' field, but it will be
        % empty.
        function key = setWithShared(storedb, store, key)
            storedb.shared = store.shared;
            store.shared = [];
            key = storedb.set(store, key);
        end
        
        % Generate a unique key and return it. This should be called
        % everytime a new point is generated / stored. Keys are valid field
        % names for structures.
        function key = getNewKey(storedb)
            key = sprintf('z%d', storedb.counter);
            storedb.counter = storedb.counter + 1;
        end
        
        % Clear entries in storedb.history to limit memory usage.
        function purge(storedb)
            
            if isinf(storedb.storedepth)
                return;
            end
            
            if storedb.storedepth <= 0
                storedb.history = struct();
                return;
            end

            % Get list of field names (keys).
            keys = fieldnames(storedb.history);
            nkeys = length(keys);

            % If we need to remove some of the elements in the database,
            if nkeys > storedb.storedepth

                % Get the last-set counter of each element:
                % a higher number means it was modified more recently.
                lastset = zeros(nkeys, 1, 'uint32');
                for i = 1 : nkeys
                    lastset(i) = storedb.history.(keys{i}).lastset__;
                end

                % Sort the counters and determine the threshold above which
                % the field needs to be removed.
                sortlastset = sort(lastset, 1, 'descend');
                minlastset = sortlastset(storedb.storedepth);

                % Remove all fields that are too old.
                storedb.history = rmfield(storedb.history, ...
                                               keys(lastset < minlastset));
            end
            
        end % end of purge()
        
    end
    
end
